# Coupled Element in ABAQUS
_Developed by Zexu Fan, Tongji University_

## 1. CPE4RUW element
_Developed based on the u-w formulated Biot equation_

Main features of this element:
- A plane-strain 2D element
- Bilinear equal-order interpolation for both the solid and fluid phases;
- Computational cost is significantly lowered using a reduced integration scheme;
- inf-sup stable (possess stability in the incompressible-impermeable limit) without additional parameters
- Interface to UMAT, thus able to directly work with existing material models

<u>It is recommended that readers cite the corresponding paper if this element (or its idea) is used in their publications, as a way to honor the authors' efforts.</u>

PS: an example concerning the dynamic consolidation of a 1-D soil column is provided in the EXP_CPE4RUW_DYNCON folder.

## 2. other new element
The 3D coupled element will be coming soon, please wait patiently.