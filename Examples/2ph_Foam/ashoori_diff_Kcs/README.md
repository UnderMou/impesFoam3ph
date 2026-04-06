# Foam Texture Sensitivity Study – Co-Injection

This case consists of four simulations with different foam texture generation rates ($K_c$), used as a base configuration inspired by the references listed below.

## How to Run

    chmod +x Allrun.sh
    ./Allrun.sh 

## Physical Description

This case considers a two-phase (water–gas) flow in porous media with foam generation during co-injection. The focus is on the impact of the **foam texture generation rate ($K_c$)** on flow behavior.

Foam is modeled through the dynamic balance between generation and destruction mechanisms, where $K_c$ controls the rate at which foam texture (i.e., lamella density) is created. Different values of $K_c$ lead to distinct foam strengths and mobility reduction effects.

Four cases are evaluated:
- Low $K_c$ → weak foam generation, higher gas mobility  
- Intermediate $K_c$ → moderate foam strength  
- High $K_c$ → strong foam generation, significant gas mobility reduction  
- Very high $K_c$ → near-local equilibrium foam behavior  

The interaction between:
- foam generation kinetics  
- phase mobility  
- viscous forces  

controls the displacement dynamics. Stronger foam reduces gas mobility, leading to improved sweep efficiency and a more stable displacement front.

This setup allows assessing the transition between **transient foam behavior** and **local equilibrium regimes**, as well as their impact on flow propagation in porous media.

## References

- Ashoori, E., Marchesin, D., & Rossen, W. R. (2010).  
  *Roles of transient and local equilibrium foam behavior in porous media – traveling wave*.  
  ECMOR XII – European Conference on the Mathematics of Oil Recovery.

- de Paula, Filipe Fernandes (2022).  
  *Conservative numerical methods to solve the two-phase flow in porous media including foam displacement*.  
  PhD thesis, Universidade Federal de Juiz de Fora, Brazil.
  <https://repositorio.ufjf.br/jspui/handle/ufjf/15152>