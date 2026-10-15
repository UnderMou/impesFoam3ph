#!/usr/bin/env python
# coding: utf-8

# # POD + non-intrusive regression ROM — `impesFoam2ph` 5-spot
# 
# Baseline: `impesFoam3ph/Examples/2ph_Well/5spot_panday`.
# 
# Only `constant/wellProperties/ConstantInjection/wells/injector1/target` is varied and interpreted as \(q_{\rm inj}\).
# 
# Pipeline:
# $
# q_{\rm inj}\to\text{FOM runs}\to S_b(\mathbf x,t;q)\to
# \text{POD}\to a(t,q)\to\text{regression}\to\hat a(t,q)\to\hat S_b.
# $
# 
# The notebook separates **POD projection error**, **regression error**, and **total ROM error**. The train/validation/test split is by complete \(q_{\rm inj}\) trajectories, avoiding leakage between times of the same simulation.
# 
# Requirements: OpenFOAM.org v9, compiled `impesFoam2ph`, `fluidfoam`, `pymor`, `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `joblib`.
# 
# > **Updated:** each sampled magnitude is written as `injector1/target = q` and `producer1/target = -q` in the copied case. The baseline case remains unchanged.
# 

# In[ ]:


# Uncomment if Python dependencies are missing:
# %pip install numpy pandas scipy matplotlib scikit-learn fluidfoam pymor joblib


# In[ ]:


from pathlib import Path
import os, re, shutil, subprocess, time, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import r2_score, mean_squared_error

from pymor.vectorarrays.numpy import NumpyVectorSpace
from pymor.algorithms.pod import pod
import fluidfoam

RANDOM_STATE = 42
rng = np.random.default_rng(RANDOM_STATE)

try:
    from IPython.display import display
except ImportError:
    display = print


# ## 1. Configuration

# In[ ]:


REPO_DIR = Path("/home/anderson/OpenFOAM/anderson-9/run/struct_impesFoam/impesFoam3ph_v6/impesFoam3ph")
print(REPO_DIR)
# REPO_URL = "https://github.com/UnderMou/impesFoam3ph.git"
BASE_CASE = REPO_DIR / "Examples/2ph_Well/5spot_panday"

SOLVER = "impesFoam2ph"
FIELD = "Sb"
FOAM_BASHRC = None  # e.g. "/opt/openfoam9/etc/bashrc" if needed

N_Q = 60
Q_FACTOR_MIN, Q_FACTOR_MAX = 0.50, 1.50
SAMPLING = "linspace"   # or "random"

CENTER_SNAPSHOTS = True
POD_ENERGY = 0.9999
MAX_POD_MODES = None

TRAIN_FRACTION = 0.67
VAL_FRACTION = 0.13

WORK_DIR = Path.cwd() / "rom_5spot_qinj"
CASES_DIR = WORK_DIR / "cases"
ARTIFACTS_DIR = WORK_DIR / "artifacts"
CASES_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

if not BASE_CASE.exists():
    print("Baseline absent; cloning public repository...")
    REPO_DIR.parent.mkdir(parents=True, exist_ok=True)
    # subprocess.run(["git", "clone", REPO_URL, str(REPO_DIR)], check=True)

assert BASE_CASE.exists(), f"Case not found: {BASE_CASE}"
print(BASE_CASE)


# ## 2. Read/change only `injector1/target`
# 
# The parser descends through `ConstantInjection → wells → injector1` and replaces only its `target`, protecting other well targets.

# In[ ]:


def find_named_block(text, name, start=0):
    m = re.search(rf"\b{re.escape(name)}\b\s*\{{", text[start:])
    if not m:
        raise ValueError(f"Block {name!r} not found")
    op = start + m.end() - 1
    depth = 0
    for i in range(op, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return op, i
    raise ValueError(f"Unclosed block {name!r}")

def well_span(text, well_name):
    """
    Locate:
    ConstantInjection
      -> wells
         -> <well_name>
    """
    lo, hi = find_named_block(text, "ConstantInjection")
    lo, hi = find_named_block(text, "wells", lo + 1)
    lo, hi = find_named_block(text, well_name, lo + 1)
    return lo, hi

def read_target(path, well_name):
    """Read target from a specific well."""
    text = Path(path).read_text()
    lo, hi = well_span(text, well_name)
    block = text[lo:hi+1]

    m = re.search(
        r"\btarget\s+([+\-0-9.eE]+)\s*;",
        block
    )
    if not m:
        raise ValueError(f"target not found for well {well_name!r}")

    return float(m.group(1))

def set_target(path, well_name, q):
    """Set and verify target for a specific well."""
    p = Path(path)
    text = p.read_text()

    lo, hi = well_span(text, well_name)
    block = text[lo:hi+1]

    block2, n = re.subn(
        r"(\btarget\s+)([+\-0-9.eE]+)(\s*;)",
        rf"\g<1>{q:.16g}\g<3>",
        block,
        count=1
    )

    if n != 1:
        raise ValueError(
            f"Could not replace target for well {well_name!r}"
        )

    p.write_text(text[:lo] + block2 + text[hi+1:])

    q_read = read_target(p, well_name)

    if not np.isclose(q_read, q, rtol=1e-12, atol=0.0):
        raise RuntimeError(
            f"Failed setting {well_name}: "
            f"requested={q}, read={q_read}"
        )

WELL_FILE = BASE_CASE / "constant/wellProperties"

q0 = read_target(WELL_FILE, "injector1")
qprod0 = read_target(WELL_FILE, "producer1")

print(f"Baseline injector target = {q0:.12e}")
print(f"Baseline producer target = {qprod0:.12e}")


# ## 3. Sample \(q_{inj}\) and split by complete simulations

# In[ ]:


if SAMPLING == "linspace":
    factors = np.linspace(Q_FACTOR_MIN, Q_FACTOR_MAX, N_Q)
else:
    factors = np.sort(rng.uniform(Q_FACTOR_MIN, Q_FACTOR_MAX, N_Q))
q_values = q0 * factors

idx = np.arange(N_Q); rng.shuffle(idx)
ntr = max(2, int(round(TRAIN_FRACTION*N_Q)))
nva = max(1, int(round(VAL_FRACTION*N_Q)))
if N_Q-ntr-nva < 1: ntr = N_Q-nva-1
split = {}
for i in idx[:ntr]: split[i] = "train"
for i in idx[ntr:ntr+nva]: split[i] = "val"
for i in idx[ntr+nva:]: split[i] = "test"

sampling_df = pd.DataFrame({
    "case_id": np.arange(N_Q), "q_inj": q_values,
    "factor_q0": factors, "split": [split[i] for i in range(N_Q)]
})
display(sampling_df.sort_values("q_inj"))


# ## 4. Generate and run OpenFOAM cases
# 
# Set `RUN_FOM=False` after the database has been generated if you only want to rebuild/retrain the ROM.

# In[ ]:


from concurrent.futures import ProcessPoolExecutor, as_completed

NP = 4
N_PARALLEL = max(1, min(NP, os.cpu_count() or 1))

print(os.cpu_count())
print(f"Running {N_Q} cases with {N_PARALLEL} parallel workers")


def is_time_dir(p):
    if not p.is_dir(): return False
    try: float(p.name); return True
    except ValueError: return False


def clean_results(case):
    for p in Path(case).iterdir():
        if is_time_dir(p) and float(p.name) > 0: shutil.rmtree(p)
        elif p.is_dir() and p.name.startswith("processor"): shutil.rmtree(p)
    for p in Path(case).glob("log.*"):
        if p.is_file(): p.unlink()


def run_shell(cmd, cwd, logfile):
    cmd = f"source {FOAM_BASHRC} && {cmd}" if FOAM_BASHRC else cmd
    with open(Path(cwd)/logfile, "w") as f:
        r = subprocess.run(["bash", "-lc", cmd], cwd=cwd,
                           stdout=f, stderr=subprocess.STDOUT, text=True)
    if r.returncode:
        raise RuntimeError(f"Failed: {cmd}; inspect {Path(cwd)/logfile}")


def prepare_case(i, q):
    dst = CASES_DIR / f"q_{i:03d}"

    if dst.exists():
        shutil.rmtree(dst)

    # Create a fresh copy. BASE_CASE is never modified.
    shutil.copytree(BASE_CASE, dst)
    clean_results(dst)

    well_file = dst / "constant/wellProperties"

    # Use the same sampled target for injector and producer.
    set_target(well_file, "injector1", q)
    set_target(well_file, "producer1", -q)

    # Explicit sanity check before running OpenFOAM.
    q_inj_read = read_target(well_file, "injector1")
    q_prod_read = read_target(well_file, "producer1")

    print(
        f"case={i:03d} | q_inj={q:.6e} | q_prod_expected={-q:.6e} | "
        f"injector={q_inj_read:.6e} | "
        f"producer={q_prod_read:.6e}"
    )

    return dst


def run_case(case):
    run_shell("blockMesh", case, "log.blockMesh")
    run_shell(SOLVER, case, f"log.{SOLVER}")


def run_one_case(case_id, q_inj, split):
    case = prepare_case(case_id, q_inj)
    print(f"[{case_id:03d}] q={q_inj:.6e} {split}")
    tic = time.perf_counter()
    run_case(case)
    dt = time.perf_counter() - tic
    print(f"  {dt:.1f} s")
    return (case_id, q_inj, split, dt)


RUN_FOM = True
records = []
if RUN_FOM:
    rows = [(int(r.case_id), float(r.q_inj), str(r.split))
            for r in sampling_df.itertuples()]

    if N_PARALLEL > 1:
        with ProcessPoolExecutor(max_workers=N_PARALLEL) as executor:
            futures = [executor.submit(run_one_case, *row) for row in rows]
            for fut in as_completed(futures):
                records.append(fut.result())
    else:
        for row in rows:
            records.append(run_one_case(*row))

    timings_df = pd.DataFrame(
        records,
        columns=["case_id", "q_inj", "split", "wall_time_s"]
    )
    display(timings_df)


# ## 5. Read common output times and `Sb` snapshots with `fluidfoam`

# In[ ]:


N_CELLS = 2601

def available_times(case):
    out={}
    for p in Path(case).iterdir():
        if is_time_dir(p): out[float(p.name)] = p.name
    return dict(sorted(out.items()))

maps=[available_times(CASES_DIR/f"q_{i:03d}") for i in range(N_Q)]
common=sorted(set.intersection(*[set(x) for x in maps]))
if not common: raise RuntimeError("No common output times")
times=np.asarray(common)
time_names=[maps[0][t] for t in common]
print(len(times), "common times:", times[:4], "...", times[-4:])

# def read_scalar(case,tname,field=FIELD):
#     return np.asarray(fluidfoam.readscalar(str(case),str(tname),field),
#                       dtype=float).squeeze().ravel()
def read_scalar(case, tname, field=FIELD):

    u = np.asarray(
        fluidfoam.readscalar(
            str(case),
            str(tname),
            field
        ),
        dtype=float
    ).squeeze().ravel()

    # OpenFOAM "internalField uniform ..."
    if u.size == 1:
        u = np.full(N_CELLS, u.item())

    if u.size != N_CELLS:
        raise ValueError(
            f"Unexpected field size in "
            f"case={case}, time={tname}, field={field}: "
            f"got {u.size}, expected {N_CELLS}"
        )

    return u

snaps=[]; rows=[]
for row in sampling_df.itertuples():
    case=CASES_DIR/f"q_{row.case_id:03d}"
    for t,tname in zip(times,time_names):
        snaps.append(read_scalar(case,tname))
        rows.append(dict(case_id=row.case_id,q_inj=row.q_inj,t=t,split=row.split))

for i, (s, row) in enumerate(zip(snaps, rows)):
    print(
        i,
        "case =", row["case_id"],
        "t =", row["t"],
        "shape =", np.asarray(s).shape
    )


X=np.vstack(snaps)
meta=pd.DataFrame(rows)
print("X:",X.shape,"(snapshots, cells)")
display(meta.groupby("split").size())


# ## 6. Finite-volume inner product
# 
# We use
# \[
# \langle u,v\rangle_V=u^TM_Vv=\sum_iV_i u_iv_i.
# $
# The notebook tries `postProcess -func writeCellVolumes`; if unavailable it falls back to equal weights, which is equivalent up to a constant on the uniform rectangular mesh.

# In[ ]:


def get_volumes(case):
    try:
        run_shell("postProcess -func writeCellVolumes -time 0",case,"log.writeCellVolumes")
        V=read_scalar(case,"0","V")
        if len(V)==X.shape[1] and np.all(V>0): return V
    except Exception as e:
        warnings.warn(str(e))
    warnings.warn("Using equal weights; appropriate up to a constant for uniform mesh.")
    return np.ones(X.shape[1])

V=get_volumes(CASES_DIR/"q_000")
sqrtV=np.sqrt(V)
print("V min/max:",V.min(),V.max(),"uniform:",np.allclose(V,V[0]))


# ## 7. POD from **training cases only**
# 
# For centered weighted snapshots
# $
# y=M_V^{1/2}(u-\bar u),
# $
# pyMOR computes Euclidean POD modes \(\Psi\). Physical modes are
# $
# \Phi=M_V^{-1/2}\Psi,\qquad \Phi^TM_V\Phi=I.
# $

# In[ ]:


train = meta.split.eq("train").to_numpy()
val   = meta.split.eq("val").to_numpy()
test  = meta.split.eq("test").to_numpy()

# -----------------------------------------
# Mean field using TRAINING snapshots only
# -----------------------------------------
u_mean = (
    X[train].mean(axis=0)
    if CENTER_SNAPSHOTS
    else np.zeros(X.shape[1])
)

# Volume-weighted centered snapshots
#
# Shape:
# Ytr = (N_train_snapshots, N_cells)
#
Ytr = (X[train] - u_mean) * sqrtV[None, :]

print("X shape    :", X.shape)
print("Ytr shape  :", Ytr.shape)
print("N cells    :", X.shape[1])
print("N train    :", train.sum())


# -----------------------------------------
# pyMOR VectorArray
# -----------------------------------------

n_cells = Ytr.shape[1]

space = NumpyVectorSpace(n_cells)

# IMPORTANT:
# this pyMOR version expects:
#
# (vector_dimension, number_of_vectors)
#
# Therefore transpose Ytr.
U = space.from_numpy(Ytr.T)

print("pyMOR VectorArray:")
print("len(U) =", len(U))
print("dim    =", U.dim)


# -----------------------------------------
# POD
# -----------------------------------------

nm = min(Ytr.shape)

if MAX_POD_MODES is not None:
    nm = min(nm, MAX_POD_MODES)

B, svals = pod(
    U,
    modes=nm
)


# -----------------------------------------
# Convert POD basis to NumPy
# -----------------------------------------

# pyMOR returns:
#
# B.to_numpy():
# (N_cells, N_modes)
#
Psi = B.to_numpy()

print("Psi shape =", Psi.shape)


# -----------------------------------------
# POD energy
# -----------------------------------------

energy = (
    np.cumsum(svals**2)
    / np.sum(svals**2)
)

r = int(
    np.searchsorted(
        energy,
        POD_ENERGY
    ) + 1
)

Psi = Psi[:, :r]


# -----------------------------------------
# Return to physical coordinates
# -----------------------------------------

Phi = Psi / sqrtV[:, None]


print()
print("Selected POD modes")
print("------------------")
print("r      =", r)
print("energy =", energy[r-1])


# -----------------------------------------
# Orthogonality check
#
# Phi^T M_V Phi = I
# -----------------------------------------

orth_error = np.linalg.norm(
    Phi.T
    @ (V[:, None] * Phi)
    - np.eye(r)
)

print(
    "orthogonality error =",
    orth_error
)


# -----------------------------------------
# Singular values
# -----------------------------------------

fig, ax = plt.subplots()

ax.semilogy(
    np.arange(1, len(svals)+1),
    svals / svals[0],
    "o-",
    ms=3
)

ax.axvline(
    r,
    ls="--"
)

ax.set(
    xlabel="POD mode",
    ylabel=r"$\sigma_i/\sigma_1$"
)

ax.grid()

plt.show()


# -----------------------------------------
# Cumulative energy
# -----------------------------------------

fig, ax = plt.subplots()

ax.plot(
    np.arange(1, len(energy)+1),
    energy
)

ax.axhline(
    POD_ENERGY,
    ls="--"
)

ax.axvline(
    r,
    ls="--"
)

ax.set(
    xlabel="Number of POD modes",
    ylabel="Cumulative energy"
)

ax.grid()

plt.show()


# ## 8. Build the target dataset $a(t,q)$
# 
# Since $\Phi^TM_V\Phi=I$,
# $
# a(t,q)=\Phi^TM_V[u(t,q)-\bar u].
# $
# In weighted coordinates this is simply `A = Y @ Psi`.

# In[ ]:


# Y=(X-u_mean)*sqrtV
A=Y@Psi
print("A:",A.shape)

coef_df=pd.concat([
    meta.reset_index(drop=True),
    pd.DataFrame(A,columns=[f"a{i+1}" for i in range(r)])
],axis=1)
display(coef_df.head())

# for j in range(min(4,r)):
#     fig,ax=plt.subplots()
#     for q in np.sort(meta.q_inj.unique()):
#         m=np.isclose(meta.q_inj,q)
#         ax.plot(meta.loc[m,"t"],A[m,j],alpha=.75)
#     ax.set(xlabel="t",ylabel=fr"$a_{j+1}$"); ax.grid(); plt.show()


import matplotlib.pyplot as plt
import matplotlib as mpl

# Unique q_inj values
q_values = np.sort(meta.q_inj.unique())

# Normalize q_inj for the colormap
norm = mpl.colors.Normalize(
    vmin=q_values.min(),
    vmax=q_values.max()
)

# Choose colormap
cmap = plt.cm.viridis

for j in range(min(4, r)):

    fig, ax = plt.subplots(figsize=(7, 5))

    for q in q_values:

        m = np.isclose(
            meta.q_inj.to_numpy(),
            q
        )

        ax.plot(
            meta.loc[m, "t"],
            A[m, j],
            color=cmap(norm(q)),
            linewidth=1.5,
            alpha=0.85
        )

    ax.set_xlabel("t")
    ax.set_ylabel(fr"$a_{{{j+1}}}$")

    ax.grid(alpha=0.3)

    # ----------------------------------
    # Colorbar representing q_inj
    # ----------------------------------

    sm = mpl.cm.ScalarMappable(
        norm=norm,
        cmap=cmap
    )

    sm.set_array([])

    cbar = fig.colorbar(
        sm,
        ax=ax
    )

    cbar.set_label(
        r"$q_{\mathrm{inj}}$"
    )

    plt.tight_layout()
    plt.show()


# ## 9. POD projection error — before any regression

# In[ ]:


def reconstruct(Ap):
    return u_mean[None,:]+Ap@Phi.T

def relL2(Utrue,Upred):
    d=np.sum(V[None,:]*(Utrue-Upred)**2,axis=1)
    n=np.sum(V[None,:]*Utrue**2,axis=1)
    return np.sqrt(d/np.maximum(n,1e-30))

Xpod=reconstruct(A)
poderr=relL2(X,Xpod)
display(pd.DataFrame({"split":meta.split,"POD_relL2":poderr})
        .groupby("split").POD_relL2.agg(["mean","median","max"]))


# ## 10. Regression: $(t,q)\to[a_1,\ldots,a_r]$
# 
# We compare Polynomial Ridge, RBF Kernel Ridge, Random Forest, Extra Trees, and MLP. The targets are standardized mode-by-mode because high-order POD coefficients may have much smaller scales.
# 
# The first implementation is multi-output. Section 13 diagnoses each \(a_i\), making it easy to switch later to one regressor per coefficient.

# In[ ]:


meta


# In[ ]:


tmin,tmax=meta.loc[train,"t"].min(),meta.loc[train,"t"].max()
def features(df):
    tau=(df.t.to_numpy()-tmin)/max(tmax-tmin,1e-30)
    eta=(df.q_inj.to_numpy()-q0)/q0
    return np.c_[tau,eta]
Z=features(meta)

Zm,Am=Z[train],A[train]
Zv,Av=Z[val],A[val]
Zt,At=Z[test],A[test]

Amean=Am.mean(0); Astd=Am.std(0); Astd[Astd<1e-14]=1
scale=lambda x:(x-Amean)/Astd
unscale=lambda x:x*Astd+Amean

models={
"Polynomial Ridge d3":make_pipeline(PolynomialFeatures(3,include_bias=False),
                                    StandardScaler(),Ridge(alpha=1e-6)),
"Kernel Ridge RBF":make_pipeline(StandardScaler(),
                                 KernelRidge(alpha=1e-6,kernel="rbf",gamma=2.0)),
"Random Forest":RandomForestRegressor(n_estimators=400,random_state=42,n_jobs=-1),
"Extra Trees":ExtraTreesRegressor(n_estimators=400,random_state=42,n_jobs=-1),
"MLP":make_pipeline(StandardScaler(),
                    MLPRegressor(hidden_layer_sizes=(64,64),activation="tanh",
                                 alpha=1e-5,max_iter=5000,early_stopping=True,
                                 random_state=42))
}


# In[ ]:


from sklearn.base import clone


class SingleOutputEnsemble:
    def __init__(self, models):
        self.models = list(models)

    def predict(self, X):
        preds = [model.predict(X) for model in self.models]
        return np.column_stack(preds)


fitted = {}
out = []

for name, model_template in models.items():
    print("training", name)
    models_per_mode = [clone(model_template) for _ in range(r)]

    for j, model in enumerate(models_per_mode):
        model.fit(Zm, scale(Am)[:, j])

    fitted[name] = SingleOutputEnsemble(models_per_mode)

    for sn, Zs, As, mask in [
        ("train", Zm, Am, train),
        ("val", Zv, Av, val),
        ("test", Zt, At, test),
    ]:
        Ahat = unscale(fitted[name].predict(Zs))
        Uh = reconstruct(Ahat)
        e = relL2(X[mask], Uh)

        out.append(dict(
            model=name,
            split=sn,
            coef_RMSE=np.sqrt(mean_squared_error(As, Ahat)),
            coef_R2=r2_score(As, Ahat, multioutput="variance_weighted"),
            field_relL2_mean=e.mean(),
            field_relL2_median=np.median(e),
            field_relL2_max=e.max(),
            Sb_min_pred=Uh.min(),
            Sb_max_pred=Uh.max(),
        ))

results = pd.DataFrame(out)
display(results.sort_values(["split", "field_relL2_mean"]))


# In[ ]:


# # MULTI-OUTPUT REGRESSOR

# fitted={}; out=[]
# for name,model in models.items():
#     print("training",name)
#     model.fit(Zm,scale(Am)); fitted[name]=model
#     for sn,Zs,As,mask in [("train",Zm,Am,train),("val",Zv,Av,val),("test",Zt,At,test)]:
#         Ah=unscale(model.predict(Zs))
#         Uh=reconstruct(Ah); e=relL2(X[mask],Uh)
#         out.append(dict(
#             model=name,split=sn,
#             coef_RMSE=np.sqrt(mean_squared_error(As,Ah)),
#             coef_R2=r2_score(As,Ah,multioutput="variance_weighted"),
#             field_relL2_mean=e.mean(),field_relL2_median=np.median(e),
#             field_relL2_max=e.max(),Sb_min_pred=Uh.min(),Sb_max_pred=Uh.max()))
# results=pd.DataFrame(out)
# display(results.sort_values(["split","field_relL2_mean"]))


# ## 11. Select by validation performance, then evaluate once on test cases

# In[ ]:


vr=results[results.split.eq("val")].sort_values("field_relL2_mean")
display(vr)
best_name=vr.iloc[0].model
best=fitted[best_name]
print("Selected:",best_name)

Ahat=unscale(best.predict(Zt))
Xhat=reconstruct(Ahat)
romerr=relL2(X[test],Xhat)

tm=meta.loc[test].reset_index(drop=True).copy()
tm["POD_projection_relL2"]=poderr[test]
tm["ROM_relL2"]=romerr
tm["regression_excess"]=romerr-poderr[test]

display(tm.groupby("q_inj")[["POD_projection_relL2","ROM_relL2","regression_excess"]]
        .agg(["mean","max"]))
print("mean POD test error:",poderr[test].mean())
print("mean ROM test error:",romerr.mean())
print("max  ROM test error:",romerr.max())


# In[ ]:


fig,ax=plt.subplots()
for q in np.sort(tm.q_inj.unique()):
    m=np.isclose(tm.q_inj,q)
    ax.plot(tm.loc[m,"t"],tm.loc[m,"ROM_relL2"],"o-",ms=3,label=f"{q:.3e}")
ax.set(xlabel="t",ylabel="ROM relative weighted L2 error")
ax.grid(); ax.legend(title="q_inj"); plt.show()


# ## 12. POD error vs ROM error
# 
# If ROM error \(\approx\) POD error, regression is already good and POD truncation dominates.  
# If ROM error \(\gg\) POD error, improve sampling/regression.  
# If both are large, global linear POD is the bottleneck.

# In[ ]:


fig,ax=plt.subplots()
ax.scatter(poderr[test],romerr)
lims=[0,max(poderr[test].max(),romerr.max())*1.05]
ax.plot(lims,lims,"--")
ax.set(xlabel="POD projection error",ylabel="total ROM error")
ax.grid(); plt.show()


# ## 13. Per-coefficient diagnostics — decide whether to use one regressor per \(a_i\)

# In[ ]:


cm=[]
for j in range(r):
    cm.append(dict(
        mode=j+1,singular_value=svals[j],cumulative_energy=energy[j],
        R2_test=r2_score(At[:,j],Ahat[:,j]) if np.var(At[:,j])>0 else np.nan,
        RMSE_test=np.sqrt(np.mean((At[:,j]-Ahat[:,j])**2))))
cm=pd.DataFrame(cm)
display(cm.head(30))

fig,ax=plt.subplots()
ax.plot(cm["mode"],cm["R2_test"],"o-")
ax.set(xlabel="POD mode",ylabel=r"test $R^2(a_i)$"); ax.grid(); plt.show()


# ## 14. Optional spatial map for the worst unseen test snapshot

# In[ ]:


iw=int(np.argmax(romerr))
utrue=X[test][iw]; urom=Xhat[iw]; ue=urom-utrue
row=tm.iloc[iw]
print(f"t={row.t}, q={row.q_inj:.6e}, ROM relL2={row.ROM_relL2:.4e}")

case=CASES_DIR/f"q_{int(row.case_id):03d}"
try:
    mesh=fluidfoam.readmesh(str(case))
    x,y=np.asarray(mesh[0]).ravel(),np.asarray(mesh[1]).ravel()
    for vals,title in [(utrue,f"FOM {FIELD}"),(urom,f"ROM {FIELD}"),(ue,"ROM - FOM")]:
        fig,ax=plt.subplots()
        sc=ax.scatter(x,y,c=vals,s=12)
        ax.set_aspect("equal"); ax.set(xlabel="x",ylabel="y",title=title)
        plt.colorbar(sc,ax=ax); plt.show()
except Exception as e:
    print("Map skipped (fluidfoam readmesh API may differ):",e)


# ## 15. Online ROM prediction and approximate speed-up

# In[ ]:


def rom_predict(t_values,q):
    df=pd.DataFrame({"t":np.asarray(t_values,dtype=float),
                     "q_inj":np.full(len(t_values),float(q))})
    Ah=unscale(best.predict(features(df)))
    return reconstruct(Ah),Ah

tic=time.perf_counter()
Unew,Anew=rom_predict(times,float(np.mean(q_values)))
romtime=time.perf_counter()-tic
print(f"{len(times)} ROM fields in {romtime:.6f} s")
if "timings_df" in globals() and len(timings_df):
    print("mean FOM run:",timings_df.wall_time_s.mean(),"s")
    print("approx speed-up:",timings_df.wall_time_s.mean()/max(romtime,1e-12))


# ## 16. Save basis, normalization, regressor, and diagnostics

# In[ ]:


np.savez_compressed(ARTIFACTS_DIR/"pod_basis.npz",
    u_mean=u_mean,Phi=Phi,V=V,singular_values=svals,energy=energy,
    q0=q0,t_min=tmin,t_max=tmax,A_mean=Amean,A_std=Astd)
joblib.dump(best,ARTIFACTS_DIR/"regressor.joblib")
sampling_df.to_csv(ARTIFACTS_DIR/"q_sampling.csv",index=False)
results.to_csv(ARTIFACTS_DIR/"model_comparison.csv",index=False)
tm.to_csv(ARTIFACTS_DIR/"test_errors.csv",index=False)
cm.to_csv(ARTIFACTS_DIR/"coefficient_metrics.csv",index=False)
print("Saved to",ARTIFACTS_DIR)


# # What to inspect first
# 
# 1. Singular-value decay and selected \(r\).
# 2. POD projection error on **test** \(q_{\rm inj}\).
# 3. Validation ranking of regressors.
# 4. Total ROM error on held-out \(q_{\rm inj}\).
# 5. \(R^2\) of each \(a_i\).
# 6. Predicted saturation range (`Sb_min_pred`, `Sb_max_pred`).
# 7. Speed-up.
# 
# Recommended next study: increase `N_Q` and plot ROM test error versus number of training \(q_{\rm inj}\) values. Then compare this multi-output approach against **one tuned regressor per POD coefficient**.
