# 2D Beam Solver

This program is designed to analyse beams by calculating their displacements under specified boundary conditions and loads. It uses FEA in combination with [LML](https://github.com/jamesbray03/Lightweight-Matrix-Library) to calculate results.

This solution is for a 2D horizontal beam, where only normal and bending loads and displacements are considered. Please click the image to see the following video on the problem.

[<img src="https://github.com/jamesbray03/FEA-with-LML/assets/47334864/fc4f82af-056e-4036-bbf5-1a9cb21d5b15">](https://youtu.be/AbDQCzjTpag?si=6DdbNPUu8S3aa9Dj)

## Usage

### Beam Properties
- **Length of Beam (meters):** This input specifies the length of the beam in meters. It defines the physical extent of the beam structure.
- **Number of Elements:** The number of elements divides the beam into smaller segments for analysis. More elements provide higher resolution but increase computational complexity.
- **Young's Modulus (GPa):** Young's modulus represents the stiffness of the material. It defines how much the material will deform under a given load.
- **Moment of Inertia (cm⁴):** Moment of inertia measures the beam's resistance to bending. Higher values indicate greater resistance to bending deformation.

### Boundary Conditions:

- **Forces:** Forces are applied along the vertical axis. Positive values represent forces acting upwards, while negative values represent forces acting downwards. Forces are measured in kilonewtons (kN).
- **Moments:** Moments represent rotational forces applied to the beam. Positive moments are anticlockwise, while negative moments are clockwise. Moments are measured in kilonewton-meters (kNm).
- **Fixpoints:** Fixpoints restrict the movement of the beam at specific nodes. Fixed supports prevent both translation and rotation, while pins allow rotation but prevent translation.

### Outputs
- **Displacements:** The program calculates and displays the vertical and rotational displacements at each node of the beam under the specified boundary conditions.
  - Vertical displacements represent the amount of vertical movement at each node of the beam in mm. 
  - Rotational displacements represent the amount of rotation at each node of the beam.

## Example

### Define the problem

![rect14099](https://github.com/jamesbray03/FEA-with-LML/assets/47334864/9a1961be-83fd-496d-b5b7-9761d43dc51e)

Let's say Young's modulus is 200GPa and the moment of intertia is 500cm4.


### Adapt for FEA

![rect14099](https://github.com/jamesbray03/FEA-with-LML/assets/47334864/ce21d814-206a-4420-b86f-84a93bccd8af)

Discreetise the beam such that conditions can be described reasonably, without a high node count. We have chosen 5 elements here.

### Use the Solver

![image](https://github.com/jamesbray03/FEA-with-LML/assets/47334864/aec59bd0-c1bb-4af3-b204-a8fb761f58a2)

![image](https://github.com/jamesbray03/FEA-with-LML/assets/47334864/7a0f3cd8-f478-4688-8c6a-8b5c32acc95b)

Here you can see the displacements and the deflections filled out for each node on the right hand side of the table.
