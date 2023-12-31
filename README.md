# Euler2D

Euler2D is a compact unstructured finite volume solver for the two-dimensional Euler equations. The code is solely based on the standard library, so no extra packages are necessary. The current version has only been tested on windows.

**Installation**

There are 2 options:

1: No installation required, simply start using the delivered .exe file

2: Compile it yourself. Since no external packages are used, installing is easy. First make sure that you have a c++ compiler installed for windows. For this project, MinGW was used in the development. The compilation is done from a makefile. So to compile it, simply open the command prompt, go to the Euler2D directory, type mingw32-make and press enter.


**Running a Simulation**

There are already 3 simulation setups given in the test directory. To run a testcase, invoke the executable followed by the location of the input file of your simulation. Example: (assuming you are in the main directory Euler2D)

' "app/euler.exe" test/shockTube/input.euler'

**How to setup a simulation**

To set up a simulation, you need to have three files in your simulation directory. You need to have an input file, a boundary file, and a mesh in vtk format. 3 examples are stored in the "test" directory of this project.

In the input file, you provide the details such as the name of the mesh and boundary file, the time step and endtime, the initial conditions, the solver settings, and the settings for saving results. 

The mesh .vtk file should contain your mesh. It can be any type of 2D cells, as long as the unstructured legacy format is used for VTK. It's recommended to use Gmsh to create your mesh, and then export it as vtk file. 

Next to the points and connectivity matrix, the mesh .vtk file should also contain cell data that assigns a value to the boundaries and the domain. This is later used to distinguish the different boundaries. A good example can be found in "test/shockTube/mesh.vtk".  

It is possible to restart simulations from a previous one, given that the mesh is exactly the same. simply use the .vtk result of your previous run as an input, see the shockTube case as an example.


In the boundary file, you associate the assigned values from the mesh.vtk file to actual physical boundaries. The options for now are sub/supersonic inlet (Dirichlet), outlet (Neumann), freestream or wall. For each boundary condition, you give the value that was assigned to the boundary mesh, and the type of boundary. 

The following variables must be declared per boundary condition:
subsonicInlet: 		u, v, rho
subsonicOutlet: 	p
supersonicInlet:	u,v,rho,p
supersonicOutlet:	None
freestream		u,v,rho,p
wall			None


**Boundary Condition Examples:**  
```
//First boundary  

Inlet  
{  
key   1			//<- key given in mesh file  
type  supersonicInlet	//<- type of boundary condition  
u     3			//<- inlet values for u,v,rho,p  
v     0  
rho   1.4  
p     1  
}  

//Second boundary  

Top Wall   
{  
key   2			//<- key given in mesh file  
type  wall		//<- type of boundary condition  
}			//no values need to be specified  

//Third boundary    
  
Outlet
{  
key   3  
type  supersonicOutlet    
}  

```


# Testcases Results

**Sod's Shocktube**  
Testing Sod's shocktube shows close results w.r.t the analytical solution. It can be noted that Euler2D shows some dissipation, which is due to the fact that it employs a first order upwind scheme, namely AUSM+.

![sods_shocktube_paraview](https://github.com/noutvdbos/Euler2D/assets/28824914/0c973cad-98a1-4c36-b6b2-e9c9b90de528)
![sods_shocktube](https://github.com/noutvdbos/Euler2D/assets/28824914/f8a5ae4b-f507-486e-849a-f2dbca791551)

**Forward Step at Mach 3**  

For the forward step case, we again see some dissipation of the shockwaves. Overall, the shockwave pattern matches closely to solutions as shown in literature. 

![forwardstep](https://github.com/noutvdbos/Euler2D/assets/28824914/7db23054-e1c5-43a3-853e-f10e961429a7)

**Naca0012 Profile at Mach 0.2**  

At Mach 0.2, a lift coefficient of 0.533 was calculated by Euler2D. According to thin airfoil theory, it should have a lift coefficient of 0.548. Giving a deviation of roughly 2 lift-counts.

![naca0012_rescaled](https://github.com/noutvdbos/Euler2D/assets/28824914/87037d1b-13aa-46cb-8cb6-71210a163e6a)
