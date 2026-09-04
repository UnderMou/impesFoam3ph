# Five-Spot Pattern – Water Injection for Oil Displacement

This case simulates a **five-spot injection pattern** for water flooding in a two-phase (water–oil) system in porous media, based on the experimental and numerical studies of Panday et al. (1994).

## How to Run

    chmod +x Allrun.sh
    ./Allrun.sh

## Physical Description

This case considers a two-phase (water–oil) flow in a porous medium using a **five-spot well pattern**, where water is injected to displace oil.

In this configuration:
- water is injected at one location (injector),  
- oil is produced at surrounding locations (producers),  
- and flow develops radially from injection toward production points.  

See figure:

<p align="center">
  <img src="schematic_5spots.pdf" alt="Five-spot pattern schematic" width="400"/>
</p>

The displacement process is governed by:
- viscous forces driving the injected water,  
- phase mobility contrast between water and oil,  
- and the geometric configuration of the well pattern.  

As water is injected, it advances through the porous medium, displacing oil toward the production wells. The efficiency of this process depends strongly on the mobility ratio and flow stability.

## Reference

- Panday, S., Wu, Y., Huyakorn, P., & Springer, E. (1994).  
  *A three-dimensional multiphase flow model for assessing NAPL contamination in porous and fractured media, 2. porous medium simulation examples*.  
  Journal of Contaminant Hydrology, 16(2), 131–156.