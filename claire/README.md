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
* **Edge_Size_First_Ortho_Wall** is the size of the first mesh edge in the first layer of blocks. This value is used to compute the refinement law (ex: = 1e-5). WARNING: this parameter will take the tiny values such as 1e-7 as equal to 0, and the algorithm will fail. You should resize your inputs.

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
Parameters not taken into account for the 2D algorithm
* **z_lim_insertions** you can let a default value. This is not use here because the algorithm works on the principle the vehicle is on plane (O,x,y), with z=0.
* **input_surface_mesh** you can let the default value (= NA), in 2D, the block discretization surface is generated automatically from the minimum number of block edges needed on the vehicle surface
* **block_surface_3D** you can let the default value (= 0)
* **nbr_linear_blocks_smoothing** you can let default value of 0, because the smoother is not available in 2D


___
### _Param.ini_ in 3D
___

Even if it is the same file as in 2D, many options are not available in 3D. As a result, you can let default parameters.
The main difference for the 3D case is that you have to give in INPUT the blocking of the surface vehicle too.

#### Dimension
* **dim** is the dimension of the mesh (=3 for the 3D case)

#### Input
* **input_mesh** is the path to the input mesh. The input mesh must be in the test_samples files. If not, you should write a path from this folder. (ex: = /Aero/2D/Apollo_2D_5.vtk). The input mesh has to be a tetrahedral mesh and the vehicle has to be completely immersed in the fluid domain.
* **input_surface_mesh** is the path to the surface blocking (ex: = /Aero/3D/C5_3D_Surface.vtk)

#### Output
* **output_file_name** is the name of the output final mesh file (ex: = AeroPipeline2D_Hex.vtk)
* **with_debug_files** this option will write intermediate files to check the state of the mesh at some steps of the algorithm (=true), or not (= false)

#### Physical parameters
* **Boundary_layer_thickness** is the thickness of the very first layer of blocks. I recommend you check carefully the size of the input domain to set a coherent size. (ex = 3.0)

#### Vector Field
* **Vector_Field_Computation_Method** choose the method to compute the vector field that will lead the direction of the extrusion (ex: = 3)
    * 0 (_default_): gradient of the combined distance field
    * 1: gradient of the distance from the vehicule
    * 2: mixt between the gradient of the distance from the vehicule (if x < **x_lim_Vector_Field_Zone1**), and the gradient of the combined distance field (if x > **x_lim_Vector_Field_Zone2**), with a damping zone in between.
* **x_lim_Vector_Field_Zone1** is the limit along the x-axis to compute the first zone of the vector field (ex: = 100)
* **x_lim_Vector_Field_Zone2** is the limit along the x-axis to compute the second zone of the vector field (ex: = 250)

#### Cell Sizes
* **Number_of_Cells_in_Boundary_layer** is the number of cells in the first layer of blocks (ex: = 30)
* **Edge_Size_on_Wall** is the default size for the edges located on the wall (ex: = 3.0)
* **Edge_Size_Default** is the default mesh size in the rest of the domain (ex: = 3.0). As we generate block-structured meshes, two opposite block edges share the same discretization. Thus, the size of some edges will be slightly larger in some part of the domain, while the size of some other edges will be smaller in other parts.
* **Edge_Size_First_Ortho_Wall** is the size of the first mesh edge in the first layer of blocks. This value is used to compute the refinement law (ex: = 1e-5). WARNING: this parameter will take the tiny values such as 1e-7 as equal to 0, and the algorithm will fail. You should resize your inputs.

#### Extrusion
* **Number_of_layers** is the number of layers (ex: =4)
* **insertions_allowed** is a boolean to say if we allow or not the use of patterns during the extrusion process.
* **insertions_allowed_on_first_layer** is a boolean to allow or not the use of patterns in the specific case of the first block layer.
* **x_lim_insertions** is the physical limit along the x-axis where you allow the insertions on a layer or not. (ex: = 100, at a node of a layer, if x<100, the insertion is not allowed, and is allowed if x>100)
* **y_lim_insertions** same as along x-axis (ex: = -10000)
* **z_lim_insertions** same as along x-axis (ex: = -10000)
* **nbr_linear_blocks_smoothing** number of specific linear block smoother iterations.

#### Curved Blocks
* **max_degree** is the maximal Bézier block degree accepted.
* **max_error** is the maximal error accepted between a block face and the vehicle surface to compute the right Bézier blocks degree. Resulting blocks can not have a Bézier block degree higher than **max_degree**, even is the condition on the error is not respected.

#### Useless parameters for the 3D
Parameters not used in the 3D algorithm
* **axisymetry**, let the default value (= false)
* **block_surface_3D** is a test/debug option used in 3D
* **Number_of_Blocks_on_Wall**, because the surface blocking is given as an input
* **Angle_of_Attack**, because corresponding vector fields not implemented yet
* **Number_of_iterations_Yao** 
* **Damping_Smoothing_Yao**
* **x_lim_SU2_inoutlet**, writing parameter for .su2 mesh file format

## /!\ WARNINGS /!\
* The number of layers can not be less than at least 2.
* **Boundary_layer_thicnkess** value has to be in the same unit as your input tri/tet mesh.
* In 2D, no matter if you want an **axisymetry** blocking as a result, and if you put the parameter to **True**, the input triangular mesh has to be of the whole fluid domain.
* Insertions/Contractions of blocks are not allowed on the **first** and **last** layer of blocks for the 2D algorithm. Therefore, if you want to generate a mesh with only two layers of blocks, there will not be any insertions.
* If you don't use the axisymmetric mode, there is no reason to have block corners on the **y axis**.
* The axisymmetric mode will cut your mesh along the **x axis**.
* For now, you will not be able to generate a final mesh with the 3D algorithm. However, you can generate a first blocking.
* Pay attention to the _param.ini_ file: it seems like if you want to set a parameter to a tiny value, it can be converted to 0 when reading in the algorithm.
Example: if _Edge_Size_First_Ortho_Wall_ is set to 1e-6, the value will be read as 1e-6. But if you decide to put 1e-7, the value read will be 0. To work with lower value, you should consider to change the scale of your input mesh.
* The methods for vector field computation 0 and 1 are inverted between the 2D and 3D algorithm.
* Allowing the use of patterns (through **insertions_allowed** for instance) does not mean patterns will be used. It still depends on the classification of the front, and the physical limits in the domain.

## ADVICES
* As a convention for the input triangular discretization of the fluid domain to mesh in 2D, the mesh is on plan (O,x,y), with z=0. We consider that the front (nose) of the vehicle is at x_min, and the back part of the vehicle is at x_max.
* As a convention for the input tetrahedral discretization of the fluid domain to mesh in 3D, we consider that the vehicle is oriented such that the vehicle front/nose is at x_min, and the vehicle back at x_max.
* The tri/tet INPUT mesh has to be refined around the geometry.
This improves the representation of the geometry in the algorithm, as we only work with a discretized representation of the geometry (and not a CAD representation of it). 
This also increases the computation of the distance field.
* In 3D, I recommand you to rescale your input geometries if they are too large (ex: if x values vary of more than 100). This impacts the accuracy of geometric classification algorithm of the input surface discretization on the geometric model given by the tetrahedral mesh. The result of a bad classification is a bad projected final hexahedral mesh.

## OUTPUT FILES 
___
#### OUTPUT 2D
___
* _AeroPipeline2D_TriMesh.vtk_ : VTK input triangular mesh, written with the fields we computed for the extrusion
* _AeroPipeline2D_QuadMesh.vtk_ : VTK file of the final mesh, written as unstructured mesh.
* _AeroPipeline2D_Blocking.vtk_ : VTK file of the linear blocking structure only.
* _AeroPipeline2D_CurvedBlocks.vtk_ : VTK trick file to visualize the curved blocking structure. This is not a real mesh. Please, use it just for visualization.
* _AeroPipeline_2D.cgns_ : Final block-structured mesh written in CGNS format

___
#### OUTPUT 3D
___
* _AeroEdgesClassification_3D\_?.vtk_ : VTK file of the block edges with the geometric classification. The ? is for the index of the front.
* _AeroExtrusion_3D\_?.vtk_ :  VTK file of the blocking already built, where i is the last set of hexa created by a pattern (regular face, or irregular edge or node). Only if **with_debug_files** is set to true in _param.ini_.
* _AeroPipeline3D_Tetra_PreTraite.vtk_ : VTK file of the input tetrahedral mesh, after being pre-treated, with the distance and vector fields computed.
* _AeroPipeline3D_Hexa.vtk_ : VTK file of the linear hex blocking.
* _Surface_3D.vtk_ : VTK file of the quad surface blocking.

## Communications linked to this work

### Paper  
* ["Block-Structured Quad Meshing for Supersonic Flow Simulations" in SIAM IMR 2023 proceedings](https://internationalmeshingroundtable.com/assets/papers/2023/11-Roche-compressed.pdf)
### Research Note  
* ["Curved Hexahedral Block-Structure Generation by Advancing Front"](https://internationalmeshingroundtable.com/assets/research-notes/imr32/2010.pdf)
### Talks  
* [Talk - SIAM IMR 2024](https://cea.hal.science/cea-04505745)
* [Talk - SIAM IMR 2023](https://hal-cea.archives-ouvertes.fr/cea-04028060)
* [Talk - 3AF 2023](https://www.3af-aerodynamics.com/images/Public/DOCS_CONFERENCE/2023/PRESENTATIONS/DAY%2001/SESSION%201B/AERO2023_38_C.%20ROCHE_PR.pdf)
### Poster
* [Poster - SIAM IMR 2023](https://hal-cea.archives-ouvertes.fr/cea-04028054)