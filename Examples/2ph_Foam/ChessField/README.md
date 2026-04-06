# Chessboard Permeability Field – Two-Phase Flow

This case simulates a two-phase (water–gas) flow in a porous medium with a **``chessboard'' permeability field**, following the benchmark configuration presented by Horgue et al. (2015).

## How to Run

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers immiscible two-phase flow in a porous medium where permeability is spatially distributed in a **checkerboard pattern**, consisting of alternating high- and low-permeability blocks.

Water and gas are injected into the domain, and the flow evolves under the combined effects of:
- strong permeability contrasts,  
- phase mobility differences,  
- and viscous forces  

The checkerboard structure creates a highly heterogeneous flow field, where fluids continuously encounter alternating regions of high and low conductivity. This leads to:

- preferential flow through high-permeability zones,  
- local flow redistribution at block interfaces,  
- distortion of saturation fronts,  
- and complex displacement patterns.  

Unlike layered systems, the heterogeneity is **multi-directional**, meaning that flow paths are not aligned with a single direction. As a result, the displacement front becomes irregular and may exhibit fingering and channeling effects depending on mobility ratios.

## Reference

- Horgue, P., Soulaine, C., Franc, J., Guibert, R., & Debenest, G. (2015).  
  *An open-source toolbox for multiphase flow in porous media*.  
  Computer Physics Communications, 187, 217–226.