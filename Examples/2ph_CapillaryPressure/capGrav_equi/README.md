# Vertical Gas–Water Column: Capillary–Gravity Equilibrium

This case simulates a vertical two-phase (water–gas) system in porous media, aiming to reach capillary and gravitational equilibrium.

## How to Run

```bash
chmod +x Allrun.sh
./Allrun.sh
```
## Physical Description

A vertical column is initially filled with water and gas. Due to density differences, gravity drives phase segregation: the lighter gas rises while the heavier water moves downward.

At the same time, capillary forces act at the interface between the phases. The competition between **gravity** and **capillary pressure** leads the system toward an equilibrium configuration.

At equilibrium:
- Gas accumulates at the top of the column  
- Water occupies the lower region  
- The pressure difference between phases satisfies the capillary pressure relation  

## Reference

Horgue, P., Soulaine, C., Franc, J., Guibert, R., & Debenest, G. (2015).  
*An open-source toolbox for multiphase flow in porous media*.  
Computer Physics Communications, 187, 217–226.

