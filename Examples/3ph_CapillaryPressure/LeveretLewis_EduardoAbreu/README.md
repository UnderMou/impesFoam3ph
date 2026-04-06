# Riemann Problems RP1 & RP2 – Capillary Diffusion Sensitivity (Leverett–Lewis Model)

This case consists of two Riemann problems (RP1 and RP2) for three-phase flow in porous media, based on the work of Abreu et al. (2004), using the Leverett–Lewis capillary pressure model with varying $\epsilon$ parameters.

## How to Run
    in each Riemman Problem (RP) directory run:

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers one-dimensional three-phase flow Riemann problems, where capillary effects are modeled using the **Leverett–Lewis capillary pressure formulation**.

Two configurations are evaluated:
- **RP1**  
- **RP2**  

Each problem represents a different initial saturation, leading to distinct wave structures during the evolution of the solution.

The main focus of this study is the impact of the parameter $\epsilon$, which controls the strength of capillary forces. Increasing $\epsilon$ increases capillary effects, introducing stronger diffusive behavior in the system.

The flow is governed by:
- viscous transport due to phase mobility  
- capillary pressure gradients  
- diffusive effects induced by capillarity  

## Reference

- Abreu, E., Furtado, F., & Pereira, F. (2004).  
  *On the Numerical Simulation of Three‐Phase Reservoir Transport Problems*.  
  Transport Theory and Statistical Physics, 33(5–7), 503–526.