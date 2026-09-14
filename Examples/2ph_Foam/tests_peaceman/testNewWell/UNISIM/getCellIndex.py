#!/usr/bin/env python3
"""
getCellIndex.py

Read an ASCII OpenFOAM polyMesh and find the cell indices associated with
one or more axis-aligned (x,y,z) well boxes.

No OpenFOAM Python bindings are required. Only NumPy is needed.

The script reconstructs the vertices belonging to each cell from:
    constant/polyMesh/points
    constant/polyMesh/faces
    constant/polyMesh/owner
    constant/polyMesh/neighbour

Two selection modes are available:

1) intersect  [default]
   Select a cell if its axis-aligned bounding box overlaps the well box.
   This is generally the most useful mode for well perforations.

2) center
   Select a cell only if the OpenFOAM cell-center approximation (mean of
   unique cell vertices) lies inside the well box.

Examples
--------
Using an external JSON file:

    python3 getCellIndex.py \
        --case . \
        --wells wells.json

Use center-based selection instead:

    python3 getCellIndex.py \
        --case . \
        --wells wells.json \
        --mode center

Write results to a JSON file:

    python3 getCellIndex.py \
        --case . \
        --wells wells.json \
        --output wellCells.json

Example wells.json
------------------
{
    "INJ1": {
        "type": "injector",
        "xmin": 1000.0,
        "xmax": 1100.0,
        "ymin": 2000.0,
        "ymax": 2100.0,
        "zmin": 2900.0,
        "zmax": 3200.0
    },

    "PROD1": {
        "type": "producer",
        "x": [5000.0, 5100.0],
        "y": [3500.0, 3600.0],
        "z": [3000.0, 3300.0]
    }
}

You can also edit DEFAULT_WELLS below and run without --wells.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Optional: define wells directly here.
# If --wells is provided, the JSON file overrides this dictionary.
# ---------------------------------------------------------------------------
DEFAULT_WELLS = {
    # "INJ1": {
    #     "xmin": 1000.0,
    #     "xmax": 1100.0,
    #     "ymin": 2000.0,
    #     "ymax": 2100.0,
    #     "zmin": 2900.0,
    #     "zmax": 3200.0,
    # },
}


def remove_comments(text: str) -> str:
    """Remove // and /* ... */ comments from an OpenFOAM ASCII file."""
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//.*?$", "", text, flags=re.M)
    return text


def extract_foam_list_body(path: Path) -> Tuple[int, str]:
    """
    Extract the main OpenFOAM list from an ASCII polyMesh file.

    Expected structure after the FoamFile dictionary:
        N
        (
            ...
        )
    """
    text = remove_comments(path.read_text())

    # Skip the FoamFile header dictionary by finding the first integer
    # followed by a parenthesized list after the header.
    matches = list(re.finditer(r"(?m)^\s*(\d+)\s*\n\s*\(", text))
    if not matches:
        raise ValueError(
            f"Could not find an OpenFOAM list in {path}. "
            "This script currently supports ASCII polyMesh files."
        )

    m = matches[-1] if len(matches) == 1 else matches[0]
    n_items = int(m.group(1))

    open_pos = text.find("(", m.start())
    if open_pos < 0:
        raise ValueError(f"Could not find opening '(' in {path}")

    # Match the top-level list parentheses.
    depth = 0
    close_pos = None
    for i in range(open_pos, len(text)):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                close_pos = i
                break

    if close_pos is None:
        raise ValueError(f"Could not find closing ')' in {path}")

    return n_items, text[open_pos + 1 : close_pos]


def read_points(path: Path) -> np.ndarray:
    """Read OpenFOAM polyMesh/points."""
    n_points, body = extract_foam_list_body(path)

    triples = re.findall(
        r"\(\s*"
        r"([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)\s*"
        r"\)",
        body,
    )

    if len(triples) != n_points:
        raise ValueError(
            f"{path}: header says {n_points} points, "
            f"but parsed {len(triples)}."
        )

    return np.asarray(triples, dtype=float)


def read_label_list(path: Path) -> np.ndarray:
    """Read owner or neighbour."""
    n_values, body = extract_foam_list_body(path)

    values = np.fromstring(body, sep=" ", dtype=np.int64)

    if len(values) != n_values:
        raise ValueError(
            f"{path}: header says {n_values} labels, "
            f"but parsed {len(values)}."
        )

    return values


def read_faces(path: Path) -> List[np.ndarray]:
    """
    Read OpenFOAM faces.

    Typical entries:
        4(0 1 2 3)
        3(4 5 6)
    """
    n_faces, body = extract_foam_list_body(path)

    face_matches = re.findall(r"(\d+)\s*\(([^()]*)\)", body)

    faces: List[np.ndarray] = []

    for n_str, labels_str in face_matches:
        expected = int(n_str)
        labels = np.fromstring(labels_str, sep=" ", dtype=np.int64)

        if len(labels) != expected:
            raise ValueError(
                f"{path}: face declares {expected} vertices "
                f"but {len(labels)} were parsed."
            )

        faces.append(labels)

    if len(faces) != n_faces:
        raise ValueError(
            f"{path}: header says {n_faces} faces, "
            f"but parsed {len(faces)}."
        )

    return faces


def reconstruct_cells(
    faces: List[np.ndarray],
    owner: np.ndarray,
    neighbour: np.ndarray,
) -> List[np.ndarray]:
    """
    Reconstruct the unique point IDs belonging to every OpenFOAM cell.

    For internal faces:
        owner[facei] and neighbour[facei]

    For boundary faces:
        only owner[facei]
    """
    if len(owner) != len(faces):
        raise ValueError(
            f"owner has {len(owner)} entries but mesh has {len(faces)} faces."
        )

    if len(neighbour) > len(faces):
        raise ValueError("neighbour list is longer than faces list.")

    max_owner = int(owner.max()) if len(owner) else -1
    max_neigh = int(neighbour.max()) if len(neighbour) else -1
    n_cells = max(max_owner, max_neigh) + 1

    cell_points = [set() for _ in range(n_cells)]

    for facei, pts in enumerate(faces):
        own = int(owner[facei])
        cell_points[own].update(map(int, pts))

        if facei < len(neighbour):
            nei = int(neighbour[facei])
            cell_points[nei].update(map(int, pts))

    return [
        np.asarray(sorted(point_ids), dtype=np.int64)
        for point_ids in cell_points
    ]


def compute_cell_geometry(
    points: np.ndarray,
    cell_points: List[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute:
        cell_min    = minimum xyz of cell vertices
        cell_max    = maximum xyz of cell vertices
        cell_center = arithmetic mean of unique cell vertices

    Note:
        cell_center is a vertex-average center, not OpenFOAM's exact
        polyhedral volume centroid. For axis-aligned box lookup it is
        normally sufficient. Use --mode intersect for perforation lookup.
    """
    n_cells = len(cell_points)

    cell_min = np.empty((n_cells, 3), dtype=float)
    cell_max = np.empty((n_cells, 3), dtype=float)
    cell_center = np.empty((n_cells, 3), dtype=float)

    for celli, pids in enumerate(cell_points):
        if len(pids) == 0:
            raise ValueError(f"Cell {celli} has no points.")

        xyz = points[pids]

        cell_min[celli] = xyz.min(axis=0)
        cell_max[celli] = xyz.max(axis=0)
        cell_center[celli] = xyz.mean(axis=0)

    return cell_min, cell_max, cell_center


def normalize_box(data: dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Accept either:

        {
            "xmin": ...,
            "xmax": ...,
            "ymin": ...,
            "ymax": ...,
            "zmin": ...,
            "zmax": ...
        }

    or:

        {
            "x": [xmin, xmax],
            "y": [ymin, ymax],
            "z": [zmin, zmax]
        }
    """
    if all(k in data for k in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")):
        lower = np.array(
            [data["xmin"], data["ymin"], data["zmin"]],
            dtype=float,
        )
        upper = np.array(
            [data["xmax"], data["ymax"], data["zmax"]],
            dtype=float,
        )

    elif all(k in data for k in ("x", "y", "z")):
        lower = np.array(
            [data["x"][0], data["y"][0], data["z"][0]],
            dtype=float,
        )
        upper = np.array(
            [data["x"][1], data["y"][1], data["z"][1]],
            dtype=float,
        )

    else:
        raise ValueError(
            "Each well must contain either "
            "xmin/xmax/ymin/ymax/zmin/zmax or x/y/z interval arrays."
        )

    # Allow the user to provide ranges in either order.
    box_min = np.minimum(lower, upper)
    box_max = np.maximum(lower, upper)

    return box_min, box_max


def find_cells_intersecting_box(
    cell_min: np.ndarray,
    cell_max: np.ndarray,
    box_min: np.ndarray,
    box_max: np.ndarray,
    tol: float = 0.0,
) -> np.ndarray:
    """
    Axis-aligned bounding-box overlap.

    A cell is selected when:
        cell_max >= box_min
    and
        cell_min <= box_max
    in x, y and z.
    """
    overlap = np.all(
        (cell_max >= box_min - tol)
        & (cell_min <= box_max + tol),
        axis=1,
    )
    return np.flatnonzero(overlap)


def find_cells_center_inside_box(
    centers: np.ndarray,
    box_min: np.ndarray,
    box_max: np.ndarray,
    tol: float = 0.0,
) -> np.ndarray:
    """Select cells whose approximate center is inside the box."""
    inside = np.all(
        (centers >= box_min - tol)
        & (centers <= box_max + tol),
        axis=1,
    )
    return np.flatnonzero(inside)


def load_wells(path: Path | None) -> Dict[str, dict]:
    if path is None:
        if not DEFAULT_WELLS:
            raise ValueError(
                "No wells supplied. Either edit DEFAULT_WELLS in the script "
                "or provide --wells wells.json."
            )
        return DEFAULT_WELLS

    with path.open() as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError("The wells JSON must contain an object/dictionary.")

    return data


def write_openfoam_label_list(path: Path, name: str, cells: np.ndarray) -> None:
    """
    Write a simple OpenFOAM label list that can be copied into dictionaries
    or adapted for well-model input.
    """
    with path.open("w") as f:
        f.write(f"// Cell indices for well {name}\n")
        f.write(f"{len(cells)}\n")
        f.write("(\n")
        for celli in cells:
            f.write(f"    {int(celli)}\n")
        f.write(")\n")


def write_well_marker_field(
    case_dir: Path,
    n_cells: int,
    well_results: Dict[str, dict],
    wells: Dict[str, dict],
    field_name: str = "wellCells",
) -> Path:
    """
    Write an OpenFOAM volScalarField:
        +1 -> injector cell
        -1 -> producer cell
         0 -> all other cells

    If two wells with conflicting signs map to the same cell, an error is raised.
    """
    values = np.zeros(n_cells, dtype=float)

    for well_name, result in well_results.items():
        well_type = str(wells[well_name].get("type", "")).strip().lower()

        if well_type in ("injector", "inj", "injection"):
            marker = 1.0
        elif well_type in ("producer", "prod", "production"):
            marker = -1.0
        else:
            raise ValueError(
                f"Well '{well_name}' has invalid or missing type. "
                "Use type='injector' or type='producer'."
            )

        ids = np.asarray(result["cell_indices"], dtype=np.int64)

        conflicts = ids[(values[ids] != 0.0) & (values[ids] != marker)]
        if len(conflicts):
            raise ValueError(
                f"Conflicting injector/producer mapping for cells "
                f"{conflicts.tolist()} involving well '{well_name}'."
            )

        values[ids] = marker

    zero_dir = case_dir / "0"
    zero_dir.mkdir(parents=True, exist_ok=True)
    out_path = zero_dir / field_name

    with out_path.open("w") as f:
        f.write(
f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  9                                     |
|   \\\\  /    A nd           |                                                 |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      {field_name};
}}

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{n_cells}
(
"""
        )

        for value in values:
            if value > 0:
                f.write("1\n")
            elif value < 0:
                f.write("-1\n")
            else:
                f.write("0\n")

        f.write(
""");

boundaryField
{
    defaultFaces
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
"""
        )

    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find OpenFOAM mesh cells associated with x/y/z well boxes."
    )

    parser.add_argument(
        "--case",
        type=Path,
        default=Path("."),
        help="OpenFOAM case directory. Default: current directory.",
    )

    parser.add_argument(
        "--wells",
        type=Path,
        default=None,
        help="JSON file containing well box coordinates.",
    )

    parser.add_argument(
        "--mode",
        choices=("intersect", "center"),
        default="intersect",
        help=(
            "Cell selection criterion. "
            "'intersect' selects cells whose bounding boxes overlap the well box; "
            "'center' selects cells whose centers lie inside it. "
            "Default: intersect."
        ),
    )

    parser.add_argument(
        "--tol",
        type=float,
        default=0.0,
        help="Coordinate tolerance added to the well box. Default: 0.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional JSON output path.",
    )

    parser.add_argument(
        "--write-lists",
        action="store_true",
        help="Also write one <wellName>_cells.txt OpenFOAM-style list per well.",
    )

    parser.add_argument(
        "--write-field",
        action="store_true",
        help=(
            "Write 0/wellCells as a volScalarField: "
            "+1 injector, -1 producer, 0 otherwise."
        ),
    )

    parser.add_argument(
        "--field-name",
        default="wellCells",
        help="Name of the marker volScalarField. Default: wellCells.",
    )

    args = parser.parse_args()

    poly_mesh = args.case / "constant" / "polyMesh"

    required = {
        "points": poly_mesh / "points",
        "faces": poly_mesh / "faces",
        "owner": poly_mesh / "owner",
        "neighbour": poly_mesh / "neighbour",
    }

    for name, path in required.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing polyMesh/{name}: {path}")

    print(f"Reading mesh: {poly_mesh}")

    points = read_points(required["points"])
    faces = read_faces(required["faces"])
    owner = read_label_list(required["owner"])
    neighbour = read_label_list(required["neighbour"])

    print(f"Points         : {len(points)}")
    print(f"Faces          : {len(faces)}")

    cells = reconstruct_cells(faces, owner, neighbour)

    print(f"Cells          : {len(cells)}")
    print("Computing cell bounding boxes...")

    cell_min, cell_max, centers = compute_cell_geometry(points, cells)

    domain_min = cell_min.min(axis=0)
    domain_max = cell_max.max(axis=0)

    print(
        "Mesh bounds    : "
        f"x=[{domain_min[0]:.8g}, {domain_max[0]:.8g}], "
        f"y=[{domain_min[1]:.8g}, {domain_max[1]:.8g}], "
        f"z=[{domain_min[2]:.8g}, {domain_max[2]:.8g}]"
    )

    wells = load_wells(args.wells)

    results = {}

    print()
    print("=" * 72)
    print(f"Selection mode: {args.mode}")
    print("=" * 72)

    for well_name, box_data in wells.items():
        box_min, box_max = normalize_box(box_data)

        if args.mode == "intersect":
            ids = find_cells_intersecting_box(
                cell_min,
                cell_max,
                box_min,
                box_max,
                tol=args.tol,
            )
        else:
            ids = find_cells_center_inside_box(
                centers,
                box_min,
                box_max,
                tol=args.tol,
            )

        results[well_name] = {
            "box_min": box_min.tolist(),
            "box_max": box_max.tolist(),
            "cell_indices": ids.astype(int).tolist(),
        }

        print()
        print(f"Well: {well_name}")
        print(
            "  box: "
            f"x=[{box_min[0]}, {box_max[0]}], "
            f"y=[{box_min[1]}, {box_max[1]}], "
            f"z=[{box_min[2]}, {box_max[2]}]"
        )
        print(f"  cells found: {len(ids)}")
        print(f"  cell indices: {ids.tolist()}")

        if args.write_lists:
            filename = Path(f"{well_name}_cells.txt")
            write_openfoam_label_list(filename, well_name, ids)
            print(f"  wrote: {filename}")

    if args.output is not None:
        with args.output.open("w") as f:
            json.dump(results, f, indent=4)

        print()
        print(f"Results written to: {args.output}")

    if args.write_field:
        field_path = write_well_marker_field(
            args.case,
            len(cells),
            results,
            wells,
            field_name=args.field_name,
        )
        print()
        print(f"Well marker field written to: {field_path}")
        print("  injector cells -> +1")
        print("  producer cells -> -1")
        print("  all other cells -> 0")


if __name__ == "__main__":
    main()
