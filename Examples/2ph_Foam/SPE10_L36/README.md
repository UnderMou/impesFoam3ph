# SPE10 – Layer 36 Full Heterogeneous Field Co-Injection

This case simulates the full heterogeneous permeability field of **Layer 36** from the Tenth SPE Comparative Solution Project in a two-phase (water–gas) co-injection scenario in porous media.

## How to Run

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers a two-phase (water–gas) flow in a porous medium characterized by a **highly heterogeneous permeability field**, corresponding to Layer 36 of the SPE10 benchmark.

Unlike simplified layered systems, the permeability field exhibits **strong spatial variability in both horizontal directions**, with permeability contrasts spanning several orders of magnitude. This complex structure represents realistic geological formations with channels, barriers, and disconnected high-permeability streaks.

Water and gas are co-injected into the domain, and the flow evolution is governed by:
- spatial permeability heterogeneity,  
- phase mobility contrast,  
- and viscous forces  

The heterogeneous field leads to:
- preferential flow through high-permeability channels,  
- flow diversion around low-permeability regions,  
- complex saturation front distortion,  
- and potential fingering and channeling effects.  

In contrast to simple layered models, flow paths are not predictable a priori and emerge from the interaction between the permeability field and multiphase flow physics.

## Reference

- Christie, M. A., & Blunt, M. J. (2001).  
  *Tenth SPE comparative solution project: A comparison of upscaling techniques*.  
  SPE Reservoir Evaluation and Engineering, 4(4), 308–317. https://doi.org/10.2118/72469-PA