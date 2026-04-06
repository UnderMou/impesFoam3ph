# Vertical Counter-Gravity Injection – 1D Two-Phase Flow

This case simulates a 1D vertical two-phase (water–gas) flow in porous media with injection in the **counter-direction of gravity**, based on the work of Durlofsky (1993).

## How to Run

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers a one-dimensional vertical flow where fluids are injected **against the direction of gravity**.

Water is injected into the domain, while gravitational forces act in the opposite direction, creating a competition between:
- viscous forces driving the injection  
- gravitational forces opposing the flow  
- phase mobility differences  

## Reference

- Durlofsky, L. J. (1993).  
  *A triangle based mixed finite element—finite volume technique for modeling two phase flow through porous media*.  
  Journal of Computational Physics, 105(2), 252–266.