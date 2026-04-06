# Secondary and Tertiary Recovery – Three-Phase Flow with Foam and Surfactant

This case consists of two simulations representing **secondary** and **tertiary recovery conditions** in a reservoir, modeled as three-phase (water–oil–gas) flow with foam and surfactant transport.

## How to Run
    In each case directory run:

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers three-phase flow in porous media including **water, oil, and gas**, with foam generation and **surfactant transport with adsorption effects**.

Two configurations are evaluated:
- **Secondary recovery case** 
- **Tertiary recovery case**  

The flow is governed by:
- multiphase Darcy flow  
- surfactant transport in the aqueous phase  
- adsorption of surfactant onto the porous matrix  
- foam generation and destruction mechanisms  

The foam model is based on an extension of the Lozano framework, where the **mobility reduction factor (MRF)** depends on surfactant concentration. In particular:
- when surfactant concentration ($C_s$) is present, **MRF > 1**, enhancing foam strength  
- when surfactant is absent, foam effects are weaker  

The interaction between:
- gas injection and foam formation  
- surfactant transport and adsorption  
- and phase mobility contrast  

controls the displacement process.

## Reference

- Mehrabi, M., Sepehrnoori, K., & Delshad, M. (2022).  
  *Displacement theory of low-tension gas flooding*.  
  Transport in Porous Media, 142(3), 475–491.

- Lozano, L. F., Chapiro, G., & Marchesin, D. (2025).  
  *The Riemann problem for three-phase foam flow in porous media*.  
  arXiv preprint arXiv:2506.08152.