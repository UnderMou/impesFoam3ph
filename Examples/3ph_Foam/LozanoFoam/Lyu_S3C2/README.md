# Three-Phase Flow with Foam – Oil Bank Formation (S3C2-Inspired Case)

This case simulates a three-phase (water–oil–gas) flow with foam in porous media, inspired by the S3C2 case of Lyu et al. (2021), using a simplified foam model based on Lozano et al. (2025).

## How to Run

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers a three-phase (water–oil–gas) flow in porous media with foam generation, aiming to reproduce **oil bank formation** during gas injection.

The configuration is inspired by a mature reservoir scenario, where gas injection (with foam effects) is used to mobilize residual oil. Foam is modeled using a simplified approach that captures the essential mechanism of **gas mobility reduction**.

The flow is governed by:
- multiphase Darcy flow (water, oil, gas)  
- foam-induced gas mobility reduction  
- phase interactions through relative permeability  

The interaction between:
- gas injection and foam formation  
- oil displacement mechanisms  
- and phase mobility contrast  

leads to the formation of an **oil bank**, where oil accumulates and is displaced ahead of the advancing gas front.

As the system evolves:
- foam reduces gas mobility, improving sweep efficiency  
- gas penetration is controlled, avoiding early breakthrough  
- oil is mobilized and concentrated into a bank  
- the oil bank propagates through the domain  

## References

- Lyu, X., Voskov, D., Tang, J., & Rossen, W. R. (2021).  
  *Simulation of foam enhanced-oil-recovery processes using operator-based linearization approach*.  
  SPE Journal, 26(04), 2287–2304.

- Lozano, L. F., Chapiro, G., & Marchesin, D. (2025).  
  *The Riemann problem for three-phase foam flow in porous media*.  
  arXiv preprint arXiv:2506.08152.