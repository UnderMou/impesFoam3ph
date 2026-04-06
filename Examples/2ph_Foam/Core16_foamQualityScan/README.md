# Foam Quality Scan – STARS Model (Valdez, 2021 Parameters)

This case reproduces a foam quality scan experiment using the STARS foam model with parameters obtained from experimental fitting performed by Valdez (2021).

## How to Run

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers a two-phase (water–gas) flow in porous media with foam generation, where the foam behavior is modeled using the **STARS implicit-texture model**.

The model parameters are taken from experimental fitting observed in core-flood experiments.

A **foam quality scan** is performed by varying the gas fraction ($f_g$) during co-injection. For each imposed foam quality, the system evolves until reaching a steady-state regime.

The simulation allows evaluating:
- the temporal evolution of pressure drop,  
- the steady-state apparent viscosity ($\mu_{app}$),  
- and the dependence of foam strength on gas fraction.  

At steady state:
- pressure drop stabilizes for each foam quality,  
- apparent viscosity can be computed from Darcy’s law,  
- and a characteristic $\mu_{app}$ vs. $f_g$ curve is obtained.  

This curve typically exhibits a peak corresponding to optimal foam strength, reflecting the balance between foam generation and collapse mechanisms.

## Outputs

This case generates:
- Pressure drop evolution over time for each foam quality  
- Steady-state apparent viscosity ($\mu_{app}$) as a function of gas fraction ($f_g$)  

## Reference

- Valdez, A. R., Rocha, B. M., da Fonseca Façanha, J. M., de Souza, A. V. O., Pérez-Gramatges, A., Chapiro, G., & Santos, R. W. D. (2022).  
  *Foam-assisted water–gas flow parameters: from core-flood experiment to uncertainty quantification and sensitivity analysis*.  
  Transport in Porous Media, 144(1), 189–209.