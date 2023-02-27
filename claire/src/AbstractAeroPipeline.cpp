//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/utils/Parameters.h>
#include <unit_test_config.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::AbstractAeroPipeline(std::string &Aparams) :
  m_meshTet(nullptr),
  m_meshHex(nullptr),
  m_manager(new cad::FACManager()),
  m_linker_TG(new cad::GeomMeshLinker()),
  m_linker_HG(new cad::GeomMeshLinker())
{
	Parameters p;
	p.add("section_Dimension","dim", Parameters::INT_P);

	p.add("section_INPUT_files","input_mesh",  Parameters::STRING_P);
	p.add("section_INPUT_files","input_surface_mesh",  Parameters::STRING_P);
	p.add("section_INPUT_files","block_surface_3D",  Parameters::INT_P);

	p.add("section_OUTPUT_files","output_file_name", Parameters::STRING_P);
	p.add("section_OUTPUT_files","with_debug_files", Parameters::BOOL_P);

	p.add("section_Physical","Boundary_layer_thickness", Parameters::DOUBLE_P);
	p.add("section_Physical","Angle_of_Attack", Parameters::DOUBLE_P);

	p.add("section_Extrusion","Number_of_layers", Parameters::INT_P);
	p.add("section_Extrusion","x_lim_insertions", Parameters::DOUBLE_P);
	p.add("section_Extrusion","y_lim_insertions", Parameters::DOUBLE_P);
	p.add("section_Extrusion","z_lim_insertions", Parameters::DOUBLE_P);

	p.add("section_Cell_Size","Number_of_Blocks_on_Wall", Parameters::INT_P);
	p.add("section_Cell_Size","Number_of_Cells_in_Boundary_Layer", Parameters::INT_P);
	p.add("section_Cell_Size","Edge_Size_on_Wall", Parameters::DOUBLE_P);
	p.add("section_Cell_Size","Edge_Size_Default", Parameters::DOUBLE_P);
	p.add("section_Cell_Size","Edge_Size_First_Ortho_Wall", Parameters::DOUBLE_P);

	p.add("section_Vector_Field","Vector_Field_Computation_Method", Parameters::INT_P);
	p.add("section_Vector_Field","x_lim_Vector_Field_Zone1", Parameters::DOUBLE_P);
	p.add("section_Vector_Field","x_lim_Vector_Field_Zone2", Parameters::DOUBLE_P);

	p.add("section_Smoothing","Number_of_iterations_Yao", Parameters::INT_P);
	p.add("section_Smoothing","Damping_Smoothing_Yao", Parameters::DOUBLE_P);

	p.add("section_SU2_Writer","x_lim_SU2_inoutlet", Parameters::DOUBLE_P);

	std::string dir(TEST_SAMPLES_DIR);
	p.parseIni(Aparams);

	p.get("section_Dimension", "dim", m_params.dimension);

	std::string input;
	p.get("section_INPUT_files","input_mesh",  input);
	m_params.input_file = dir+input;
	p.get("section_INPUT_files","input_surface_mesh",  input);
	m_params.input_file_3D_surface = dir+input;
	p.get("section_INPUT_files","block_surface_3D",  m_params.block_surface_3D);

	p.get("section_OUTPUT_files","output_file_name", m_params.output_file);
	p.get("section_OUTPUT_files","with_debug_files", m_params.with_debug_files);

	p.get("section_Physical","Boundary_layer_thickness", m_params.delta_cl);
	p.get("section_Physical","Angle_of_Attack", m_params.angle_attack);

	p.get("section_Extrusion","Number_of_layers", m_params.nbr_couches);
	p.get("section_Extrusion","x_lim_insertions", m_params.x_lim);
	p.get("section_Extrusion","y_lim_insertions", m_params.y_lim);
	p.get("section_Extrusion","z_lim_insertions", m_params.z_lim);

	p.get("section_Cell_Size","Number_of_Blocks_on_Wall", m_params.nbrMinBloc);
	p.get("section_Cell_Size","Number_of_Cells_in_Boundary_Layer", m_params.nbrCellsInCL);
	p.get("section_Cell_Size","Edge_Size_on_Wall", m_params.edge_size_wall);
	p.get("section_Cell_Size","Edge_Size_Default", m_params.edge_size_default);
	p.get("section_Cell_Size","Edge_Size_First_Ortho_Wall", m_params.edge_size_first_ortho_wall);

	p.get("section_Vector_Field","Vector_Field_Computation_Method", m_params.vectors_field);
	p.get("section_Vector_Field","x_lim_Vector_Field_Zone1", m_params.x_VectorField_Z1);
	p.get("section_Vector_Field","x_lim_Vector_Field_Zone2", m_params.x_VectorField_Z2);

	p.get("section_Smoothing","Number_of_iterations_Yao", m_params.nbr_iter_smoothing_yao);
	p.get("section_Smoothing","Damping_Smoothing_Yao", m_params.damping_smoothing_yao);

	p.get("section_SU2_Writer","x_lim_SU2_inoutlet", m_params.x_lim_SU2_inoutlet);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::~AbstractAeroPipeline()
{
	delete m_manager;
	delete m_linker_TG;
	delete m_linker_HG;
	delete m_meshTet;
	delete m_meshHex;
}
/*------------------------------------------------------------------------*/