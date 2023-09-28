# Aero module

This module aims to provide a solution to generate meshes dedicated to Computational Fluid Dynamics.
You can find more details about the 2D method in [this paper](https://internationalmeshingroundtable.com/assets/papers/2023/11-Roche-compressed.pdf).

## How to run a test with the algorithm

You can run the full pipeline via the file _AeroPipelineTestSuite.h_ in the tst folder, for both 2D and 3D cases.
To set the parameters for the algorithm, you can fulfill a file _param.ini_.

___
### _Param.ini_ in 2D
___
I recommend you check the file _param_Apollo_2D.ini_.
All the input files need to be VTK files.

#### Dimension
* **dim** is the dimension of the mesh (=2 for the 2D case)
* **axisymetry** is a boolean to choose if the output mesh should be axis-symmetric (= true) or not (= false)  
  
#### Input
* **input_mesh** is the path to the input mesh. The input mesh must be in the test_samples files. If not, you should write a path from this folder. (ex: = /Aero/2D/Apollo_2D_5.vtk). The input mesh has to be a triangular mesh and the vehicle has to be completely immersed in the fluid surface. You can not give as an input a triangular mesh on an axis-symmetric domain.

#### Output
* **output_file_name** is the name of the output final mesh file (ex: = AeroPipeline2D_Quad.vtk)
* **with_debug_files** this option will write intermediate files to check the state of the mesh at some steps of the algorithm (=true), or not (= false)

#### Physical parameters
* **Boundary_layer_thickness** is the thickness of the very first layer of blocks. I recommend you check carefully the size of the input domain to set a coherent size. (ex = 3.0)
* **Angle_of_Attack** is an angle in degrees used to compute the vector field for the extrusion. (ex: = 0.0)

#### Cell Sizes
* **Number_of_Blocks_on_Wall** is the minimum number of block edges on the wall (ex: = 4)
* **Number_of_Cells_in_Boundary_layer** is the number of cells in the first layer of blocks (ex: = 30)
* **Edge_Size_on_Wall** is the default size for the edges located on the wall (ex: = 3.0)
* **Edge_Size_Default** is the default mesh size in the rest of the domain (ex: = 3.0). As we generate block-structured meshes, two opposite block edges share the same discretization. Thus, the size of some edges will be slightly larger in some part of the domain, while the size of some other edges will be smaller in other parts.
* **Edge_Size_First_Ortho_Wall** is the size of the first mesh edge in the first layer of blocks. This value is used to compute the refinement law (ex: = 1e-5). ATTENTION: this parameter will take the tiny values such as 1e-7 as equal to 0, and the algorithm will fail. You should resize your inputs.

#### Extrusion
* **Number_of_layers** is the number of layers (ex: =4)
* **x_lim_insertions** is the physical limit along the x-axis where you allow the insertions on a layer or not. (ex: = 100, at a node of a layer, if x<100, the insertion is not allowed, and is allowed if x>100)
* **y_lim_insertions** same as along x-axis (ex: = -10000)
* **z_lim_insertions** same as along x-axis (ex: = -10000)

The insertions are not allowed on the last layer of blocks, to avoid flattened blocks on the far-field boundary. So if you choose to have only two layers of blocks, there will be no insertions.

#### Vector Field

* **Vector_Field_Computation_Method** choose the method to compute the vector field that will lead the direction of the extrusion (ex: = 3)
  * 0 (_default_): gradient of the distance from the vehicule
  * 1: gradient of the combined distance field
  * 2: mixt between the gradient of the distance from the vehicule (if x < **x_lim_Vector_Field_Zone1**), and the gradient of the combined distance field (if x > **x_lim_Vector_Field_Zone2**), with a damping zone in between.
  * 3: mixt between the gradient of the distance from the vehicule (if x < **x_lim_Vector_Field_Zone1**), and the angle of attack (if x > **x_lim_Vector_Field_Zone2**), with a damping zone in between.
  * 4: mixt between the gradient of the combined distance field (if x < **x_lim_Vector_Field_Zone1**), and the angle of attack (if x > **x_lim_Vector_Field_Zone2**), with a damping zone in between.
* **x_lim_Vector_Field_Zone1** is the limit along the x-axis to compute the first zone of the vector field (ex: = 100)
* **x_lim_Vector_Field_Zone2** is the limit along the x-axis to compute the second zone of the vector field (ex: = 250)

I highly recommend you to check on the [IMR paper](https://internationalmeshingroundtable.com/assets/papers/2023/11-Roche-compressed.pdf) to understand well this section.
You'll find some pictures to understand the differents parameters.

#### Smoothing
These are the parameters for the Yao Smoothing on the first layer of blocks.

* **Number_of_iterations_Yao** is the number of iterations for the Yao Smoother (ex = 10)
* **Damping_Smoothing_Yao** is the damping parameter, between 0 and 1 (ex = 0.2). If 0, no damping, if 1, the mesh will not be smoothed.

#### Useless parameters for the 2D
* **input_surface_mesh** you can let the default value (= NA)
* **block_surface_3D** you can let the default value (= 0)


___
### _Param.ini_ in 3D
___

Even if it is the same file as in 2D, many options are not available in 3D. As a result, you can let default parameters.
The main difference for the 3D case is that you have to give in INPUT the blocking of the surface vehicle too.

#### Surface blocking
* **input_surface_mesh** is the path to the surface blocking (ex: = /Aero/3D/C5_3D_Surface.vtk)
* **block_surface_3D**, let equal to the default value 0 if you give the surface blocking as an INPUT.

#### Vector Field
* **Vector_Field_Computation_Method** choose the method to compute the vector field that will lead the direction of the extrusion (ex: = 3)
    * 0 (_default_): gradient of the combined distance field
    * 1: gradient of the distance from the vehicule
    * 2: mixt between the gradient of the distance from the vehicule (if x < **x_lim_Vector_Field_Zone1**), and the gradient of the combined distance field (if x > **x_lim_Vector_Field_Zone2**), with a damping zone in between.
* **x_lim_Vector_Field_Zone1** is the limit along the x-axis to compute the first zone of the vector field (ex: = 100)
* **x_lim_Vector_Field_Zone2** is the limit along the x-axis to compute the second zone of the vector field (ex: = 250)

#### Useless parameters for the 3D
* **axisymetry**, let the default value (= false)
* **Number_of_Blocks_on_Wall**, because the surface blocking is given as an input
* **x_lim_insertions**
* **y_lim_insertions**
* **z_lim_insertions**
* **Number_of_iterations_Yao** 
* **Damping_Smoothing_Yao**

### INPUT Tri/Tet mesh

## WARNINGS
* The number of layers can not be less than at least 2.
* Insertions/Contractions of blocks are not allowed on the **first** and **last** layer of blocks for the 2D algorithm. Therefore, if you want to generate a mesh with only two layers of blocks, there will not be any insertions.
* If you don't use the axisymmetric mode, there is no reason to have block corners on the **y axis**.
* For now, you will not be able to generate a final mesh with the 3D algorithm. However, you can generate a first blocking.
* Pay attention to the _param.ini_ file: it seems like if you want to set a parameter to a tiny value, it can be converted to 0 when reading in the algorithm.
Example: if _Edge_Size_First_Ortho_Wall_ is set to 1e-6, the value will be read as 1e-6. But if you decide to put 1e-7, the value read will be 0.
* The methods for vector field computation 0 and 1 are inverted between the 2D and 3D algorithm.

## Advices
* The tri/tet INPUT mesh has to be refined around the geometry.
This improves the representation of the geometry in the algorithm, as we only work with a discretized representation of the geometry. 
This also increases the computation of the distance field.

## OUTPUT
* _AeroPipeline2D_QuadMesh.vtk_ : .vtk file of the final mesh, written as unstructured mesh.
* _AeroPipeline2D_Blocking.vtk_ : .vtk file of the linear blocking structure only.
* _AeroPipeline2D_CurvedBlocks.vtk_ : .vtk file of the curved blocking structure.

## Communications linked to this work

* ["Block-Structured Quad Meshing for Supersonic Flow Simulations" in SIAM IMR 2023 proceedings](https://internationalmeshingroundtable.com/assets/papers/2023/11-Roche-compressed.pdf)
* [Talk - SIAM IMR 2023](https://hal-cea.archives-ouvertes.fr/cea-04028060)
* [Poster - SIAM IMR 2023](https://hal-cea.archives-ouvertes.fr/cea-04028054)
* [Talk - 3AF 2023](https://www.3af-aerodynamics.com/images/Public/DOCS_CONFERENCE/2023/PRESENTATIONS/DAY%2001/SESSION%201B/AERO2023_38_C.%20ROCHE_PR.pdf)