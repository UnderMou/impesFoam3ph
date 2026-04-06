# Gas–Water co-injection with Foam and Surfactant Transport

This test case contains two cases simulating a co-injection two-phase (water–gas) system in porous media using the **STARS foam model** with surfactant transport and adsorption.

## How to Run
    in each case:

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

Each directory in this repository corresponds to a simulation of a gas–water co-injection, where foam is generated during flow and its behavior is influenced by **surfactant transport and adsorption**.

The physical system also accounts for:  
- surfactant transport in the aqueous phase,  
- and adsorption of surfactant onto the solid matrix.  

Surfactant plays a key role in stabilizing foam by reducing gas mobility. Its transport and adsorption directly affect local foam strength, introducing additional coupling between flow and chemical processes.

## Reference

- Fritis, G. C., Paz, P. S., Lozano, L. F., & Chapiro, G. (2024).  
  *On the Riemann problem for the foam displacement in porous media with linear adsorption*.  
  SIAM Journal on Applied Mathematics, 84(2), 581–601.

