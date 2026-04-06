# Two-layeres heterogeneity permeability field co-injection

This case simulates a 2-Layered (different K between layers) in a two-phase (water–gas) system in porous media co-injection.

## How to Run

```bash
chmod +x Allrun.sh
./Allrun.sh
```
## Physical Description

This case considers a two-phase (water–gas) flow in a porous medium composed of two horizontal layers with different permeabilities.

Water and gas are co-injected into the domain, and their flow is strongly influenced by the permeability contrast between the layers. The higher-permeability layer offers less resistance, promoting preferential flow and faster phase transport, while the lower-permeability layer restricts flow.

The interaction between:
- permeability heterogeneity,  
- phase mobility,  
- and viscous forces  

leads to an uneven distribution of saturation, with channeling effects typically occurring in the high-permeability layer.

Depending on the relative permeabilities and injection conditions, crossflow between layers may also occur, redistributing fluids and affecting the overall displacement efficiency.

This setup is representative of layered reservoirs. 

## Reference

no reference

