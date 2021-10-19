/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Parameters.cpp
 *  \author  legoff
 *  \date    04/24/2019
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/Parameters.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Project File Headers
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
int Parameters::dimension = -1;
/*----------------------------------------------------------------------------*/
kmds::ECellType Parameters::celltype = kmds::KMDS_UNDEF;
/*----------------------------------------------------------------------------*/
bool Parameters::badpillowing_detection_ON = true;
/*----------------------------------------------------------------------------*/
Parameter_badpillowing_method Parameters::badpillowing_detection_method = Parameter_badpillowing_method::glpk;
/*----------------------------------------------------------------------------*/
int Parameters::pillow_node_seed_range = 2;
/*----------------------------------------------------------------------------*/
double Parameters::quality_threshold = 0.3;
/*----------------------------------------------------------------------------*/
std::string  Parameters::output_prefix = "./";
/*----------------------------------------------------------------------------*/
int Parameters::debug_verbosity_lvl = 1;
/*----------------------------------------------------------------------------*/
bool Parameters::debug_ouputfile_pillow = false;
/*----------------------------------------------------------------------------*/
bool Parameters::debug_ouputfile_pixels = true;
/*----------------------------------------------------------------------------*/
std::map<std::pair<int, int>, int> Parameters::mat_pair_2_int = {{std::pair<int, int>(0,1),  0},
                                                                 {std::pair<int, int>(0,2),  1},
                                                                 {std::pair<int, int>(1,2),  2},
                                                                 {std::pair<int, int>(0,3),  3},
                                                                 {std::pair<int, int>(1,3),  4},
                                                                 {std::pair<int, int>(2,3),  5},
    };

/*----------------------------------------------------------------------------*/
// 3-refinement
std::set<tableMarkedNodes<4>, compare_tableMarkedNodes<4> > Refinement_tables::lookupNodesValid_3_2D =
    {
            { tableMarkedNodes<4> { 0,0,0,0} }
    };

std::map<tableMarkedNodes<4>, tableMarkedNodes<4>, compare_tableMarkedNodes<4> > Refinement_tables::lookupNodesSolve_3_2D =
    {
            { std::pair<tableMarkedNodes<4>, tableMarkedNodes<4> > ( tableMarkedNodes<4> { 0,0,0,0}, tableMarkedNodes<4> { 0, 0, 0, 0})}
    };

std::set<tableMarkedNodes<8>, compare_tableMarkedNodes<8> > Refinement_tables::lookupNodesValid_3_3D =
        {
                // all
                { tableMarkedNodes<8> { 0,0,0,0,0,0,0,0} },
                { tableMarkedNodes<8> { 1,1,1,1,1,1,1,1} },

                // corners
                { tableMarkedNodes<8> { 1,0,0,0,0,0,0,0} },
                { tableMarkedNodes<8> { 0,1,0,0,0,0,0,0} },
                { tableMarkedNodes<8> { 0,0,1,0,0,0,0,0} },
                { tableMarkedNodes<8> { 0,0,0,1,0,0,0,0} },
                { tableMarkedNodes<8> { 0,0,0,0,1,0,0,0} },
                { tableMarkedNodes<8> { 0,0,0,0,0,1,0,0} },
                { tableMarkedNodes<8> { 0,0,0,0,0,0,1,0} },
                { tableMarkedNodes<8> { 0,0,0,0,0,0,0,1} },

                // edges
                { tableMarkedNodes<8> { 1,1,0,0,0,0,0,0} },
                { tableMarkedNodes<8> { 0,1,1,0,0,0,0,0} },
                { tableMarkedNodes<8> { 0,0,1,1,0,0,0,0} },
                { tableMarkedNodes<8> { 1,0,0,1,0,0,0,0} },
                { tableMarkedNodes<8> { 1,0,0,0,1,0,0,0} },
                { tableMarkedNodes<8> { 0,1,0,0,0,1,0,0} },
                { tableMarkedNodes<8> { 0,0,1,0,0,0,1,0} },
                { tableMarkedNodes<8> { 0,0,0,1,0,0,0,1} },
                { tableMarkedNodes<8> { 0,0,0,0,1,1,0,0} },
                { tableMarkedNodes<8> { 0,0,0,0,0,1,1,0} },
                { tableMarkedNodes<8> { 0,0,0,0,0,0,1,1} },
                { tableMarkedNodes<8> { 0,0,0,0,1,0,0,1} },

                // faces
                { tableMarkedNodes<8> { 1,1,1,1,0,0,0,0} },
                { tableMarkedNodes<8> { 0,0,0,0,1,1,1,1} },
                { tableMarkedNodes<8> { 1,1,0,0,1,1,0,0} },
                { tableMarkedNodes<8> { 0,0,1,1,0,0,1,1} },
                { tableMarkedNodes<8> { 1,0,0,1,1,0,0,1} },
                { tableMarkedNodes<8> { 0,1,1,0,0,1,1,0} },
        };

std::map<tableMarkedNodes<8>, tableMarkedNodes<8>, compare_tableMarkedNodes<8> > Refinement_tables::lookupNodesSolve_3_3D =
        {
                // all
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,0,0,0}, tableMarkedNodes<8> { 0,0,0,0,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,1,1,1,1,1,1}, tableMarkedNodes<8> { 1,1,1,1,1,1,1,1}) },

                // corners
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,0,0,0,0,0}, tableMarkedNodes<8> { 1,0,0,0,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,0,0,0,0,0,0}, tableMarkedNodes<8> { 0,1,0,0,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,0,0,0,0,0}, tableMarkedNodes<8> { 0,0,1,0,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,1,0,0,0,0}, tableMarkedNodes<8> { 0,0,0,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,0,0,0}, tableMarkedNodes<8> { 0,0,0,0,1,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,1,0,0}, tableMarkedNodes<8> { 0,0,0,0,0,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,0,1,0}, tableMarkedNodes<8> { 0,0,0,0,0,0,1,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,0,0,1}, tableMarkedNodes<8> { 0,0,0,0,0,0,0,1}) },

                // edges
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,0,0,0,0,0,0}, tableMarkedNodes<8> { 1,1,0,0,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,1,0,0,0,0,0}, tableMarkedNodes<8> { 0,1,1,0,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,1,0,0,0,0}, tableMarkedNodes<8> { 0,0,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,1,0,0,0,0}, tableMarkedNodes<8> { 1,0,0,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,0,1,0,0,0}, tableMarkedNodes<8> { 1,0,0,0,1,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,0,0,0,1,0,0}, tableMarkedNodes<8> { 0,1,0,0,0,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,0,0,0,1,0}, tableMarkedNodes<8> { 0,0,1,0,0,0,1,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,1,0,0,0,1}, tableMarkedNodes<8> { 0,0,0,1,0,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,1,0,0}, tableMarkedNodes<8> { 0,0,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,1,1,0}, tableMarkedNodes<8> { 0,0,0,0,0,1,1,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,0,1,1}, tableMarkedNodes<8> { 0,0,0,0,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,0,0,1}, tableMarkedNodes<8> { 0,0,0,0,1,0,0,1}) },

                // diagonals
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,1,0,0,0,0,0}, tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,0,1,0,0,0,0}, tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,0,1,0}, tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,1,0,1}, tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,0,0,1,0,0}, tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,0,0,1,0,0,0}, tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,0,0,0,0,1}, tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,1,0,0,1,0}, tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,1,1,0,0,0}, tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,0,0,0,0,1}, tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,0,0,1,0,0}, tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,0,0,0,0,1,0}, tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}) },

                // faces
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}, tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}, tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}, tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}, tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}, tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}, tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}) },

                // faces 3 nodes
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,1,1,0,0,0,0}, tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,1,1,0,0,0,0}, tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,0,1,0,0,0,0}, tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,1,0,0,0,0,0}, tableMarkedNodes<8> { 1,1,1,1,0,0,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,0,1,1,1}, tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,0,1,1}, tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,1,0,1}, tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,0,1,1,1,0}, tableMarkedNodes<8> { 0,0,0,0,1,1,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,0,0,1,1,0,0}, tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,0,1,1,0,0}, tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,0,0,0,1,0,0}, tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,1,0,0,1,0,0,0}, tableMarkedNodes<8> { 1,1,0,0,1,1,0,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,1,0,0,1,1}, tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,0,0,0,1,1}, tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,1,0,0,0,1}, tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,1,0,0,1,0}, tableMarkedNodes<8> { 0,0,1,1,0,0,1,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,0,1,1,0,0,1}, tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,0,1,0,0,1}, tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,1,0,0,0,1}, tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 1,0,0,1,1,0,0,0}, tableMarkedNodes<8> { 1,0,0,1,1,0,0,1}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,0,1,0,0,1,1,0}, tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,0,0,0,1,1,0}, tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,1,0,0,0,1,0}, tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}) },
                { std::pair<tableMarkedNodes<8>, tableMarkedNodes<8> > ( tableMarkedNodes<8> { 0,1,1,0,0,1,0,0}, tableMarkedNodes<8> { 0,1,1,0,0,1,1,0}) },
        };

    std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > Refinement_tables::nbNodesOnFaceToBuild_3a_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,0},  0) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,1},  4) },

                    // corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,1},  1) },

                    // edges
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,1},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,1},  4) },

                    // diagonals
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,1},  1) },

                    // three-corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,1},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,1},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,1},  4) },
            };

    std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > Refinement_tables::nbNodesOnFaceToBuild_3b_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,0},  0) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,1},  4) },

                    // corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,1},  1) },

                    // edges
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,1},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,1},  4) },

                    // diagonals
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,1},  4) },

                    // three-corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,1},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,1},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,1},  4) },
            };

    std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > Refinement_tables::nbNodesOnFaceToBuild_3b_3D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,0},  0) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,1}, 12) },

                    // corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,1},  3) },

                    // edges
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,0},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,0},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,1},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,1},  8) },
            };

    std::array<int, 8> Refinement_tables::nodesOnFaceFromExternal_3ab_2D =
            {
                    2*7 + 0,  // bot
                    4*7 + 0,
                    6*7 + 2,  // right
                    6*7 + 4,
                    2*7 + 6,  // top
                    4*7 + 6,
                    0*7 + 2,  // left
                    0*7 + 4
            };

    std::map<tableMarkedNodes<4>, std::vector<std::array<int, 2> >, compare_tableMarkedNodes<4> > Refinement_tables::nodesOnFaceToBuild_3a_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 0, 0},
                            {

                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},

                    // corners
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 0, 0},
                            {
                                    {2,2},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 0, 0},
                            {
                                    {4,2},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 1, 0},
                            {
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 0, 1},
                            {
                                    {2,4},
                            })},

                    // edges
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 0, 0},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 1, 0},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 0, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},

                    // diagonals
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 1, 0},
                            {
                                    {3,3},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 0, 1},
                            {
                                    {3,3},
                            })},

                    // three-corners
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 1, 0},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 0, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
            };

    std::map<tableMarkedNodes<4>, std::vector<std::array<int, 2> >, compare_tableMarkedNodes<4> > Refinement_tables::nodesOnFaceToBuild_3b_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 0, 0},
                            {

                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},

                    // corners
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 0, 0},
                            {
                                    {2,2},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 0, 0},
                            {
                                    {4,2},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 1, 0},
                            {
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 0, 1},
                            {
                                    {2,4},
                            })},

                    // edges
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 0, 0},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 1, 0},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 0, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 0, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},

                    // diagonals
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 1, 0},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 0, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},

                    // three-corners
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 1, 0},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {0, 1, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 0, 1, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
                    { std::pair<tableMarkedNodes<4>, std::vector<std::array<int, 2> > > (
                            tableMarkedNodes<4> {1, 1, 0, 1},
                            {
                                    {2,2},
                                    {4,2},
                                    {2,4},
                                    {4,4},
                            })},
            };

    std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > Refinement_tables::nbFacesOnFaceToBuild_3a_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,1},  9) },

                    // corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,1},  3) },

                    // edges
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,0},  7) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,0},  7) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,1},  7) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,1},  7) },

                    // diagonals
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,0},  4) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,1},  4) },

                    // three-corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,0},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,1},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,1},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,1},  8) },
            };

    std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > Refinement_tables::nbFacesOnFaceToBuild_3b_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,0},  1) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,1},  9) },

                    // corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,0},  3) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,0,1},  3) },

                    // edges
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,0},  7) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,0},  7) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,0,1,1},  7) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,0,1},  7) },

                    // diagonals
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,0},  7) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,0,1},  7) },

                    // three-corners
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,1,0},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {0,1,1,1},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,0,1,1},  8) },
                    { std::pair<tableMarkedNodes<4>, int > ( tableMarkedNodes<4> {1,1,0,1},  8) },
            };

    std::map<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > Refinement_tables::facesOnFaceToBuild_3a_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,0,0},
                                                                                           {
                                                                                                   { 0, 6,48,42}
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,1,1},
                                                                                           {
                                                                                                   { 0, 2,16,14},
                                                                                                   { 2, 4,18,16},
                                                                                                   { 4, 6,20,18},
                                                                                                   {14,16,30,28},
                                                                                                   {16,18,32,30},
                                                                                                   {18,20,34,32},
                                                                                                   {28,30,44,42},
                                                                                                   {30,32,46,44},
                                                                                                   {32,34,48,46},
                                                                                           } ) },

                    // corners
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,0,0},
                                                                                           {
                                                                                                   { 0, 2,16,14},
                                                                                                   { 2, 6,48,16},
                                                                                                   {14,16,48,42},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,0,0},
                                                                                           {
                                                                                                   { 6,20,18, 4},
                                                                                                   {20,48,42,18},
                                                                                                   { 4,18,42, 0},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,1,0},
                                                                                           {
                                                                                                   {48,46,32,34},
                                                                                                   {46,42, 0,32},
                                                                                                   {34,32, 0, 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,0,1},
                                                                                           {
                                                                                                   {42,28,30,44},
                                                                                                   {28, 0, 6,30},
                                                                                                   {44,30, 6,48},
                                                                                           } ) },

                    // edges
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,0,0},
                                                                                           {
                                                                                                   { 0, 2,16,14},
                                                                                                   { 2, 4,18,16},
                                                                                                   { 4, 6,20,18},
                                                                                                   {14,16,30,42},
                                                                                                   {16,18,32,30},
                                                                                                   {18,20,48,32},
                                                                                                   {30,32,48,42},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,1,0},
                                                                                           {
                                                                                                   { 6,20,18, 4},
                                                                                                   {20,34,32,18},
                                                                                                   {34,48,46,32},
                                                                                                   { 4,18,16, 0},
                                                                                                   {18,32,30,16},
                                                                                                   {32,46,42,30},
                                                                                                   {16,30,42, 0},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,1,1},
                                                                                           {
                                                                                                   {48,46,32,34},
                                                                                                   {46,44,30,32},
                                                                                                   {44,42,28,30},
                                                                                                   {34,32,18, 6},
                                                                                                   {32,30,16,18},
                                                                                                   {30,28, 0,16},
                                                                                                   {18,16, 0, 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,0,1},
                                                                                           {
                                                                                                   {48,46,32,34},
                                                                                                   {46,44,30,32},
                                                                                                   {44,42,28,30},
                                                                                                   {34,32,18, 6},
                                                                                                   {32,30,16,18},
                                                                                                   {30,28, 0,16},
                                                                                                   {18,16, 0, 6},
                                                                                           } ) },
//                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,0,1},  7) },
//
//                    // diagonals
//                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,1,0},  4) },
//                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,0,1},  4) },
//
//                    // three-corners
//                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,1,0},  8) },
//                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,1,1},  8) },
//                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,1,1},  8) },
//                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,0,1},  8) },
            };

    std::map<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > Refinement_tables::facesOnFaceToBuild_3b_2D =
            {
                    // all
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,0,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 6*7 + 0, 6*7 + 6, 0*7 + 6}
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,1,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 4*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 2*7 + 4, 0*7 + 4},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 2, 6*7 + 4, 4*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 4*7 + 4, 4*7 + 6, 2*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },

                    // corners
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,0,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 6*7 + 0, 6*7 + 6, 2*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 6*7 + 6, 0*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,0,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 4*7 + 0, 4*7 + 2, 0*7 + 6},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 6, 4*7 + 2, 6*7 + 2, 6*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,1,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 6*7 + 0, 6*7 + 4, 4*7 + 4},
                                                                                                   { 0*7 + 0, 4*7 + 4, 4*7 + 6, 0*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,0,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 6*7 + 0, 2*7 + 4, 0*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 6*7 + 0, 6*7 + 6, 2*7 + 6},
                                                                                           } ) },

                    // edges
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,0,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 4*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 2*7 + 4, 0*7 + 6},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 2, 6*7 + 6, 4*7 + 4},
                                                                                                   { 2*7 + 4, 4*7 + 4, 6*7 + 6, 0*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,1,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 4*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 0, 2*7 + 2, 2*7 + 4, 0*7 + 6},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 2, 6*7 + 4, 4*7 + 4},
                                                                                                   { 2*7 + 4, 4*7 + 4, 4*7 + 6, 0*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,0,1,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 6*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 0*7 + 0, 2*7 + 2, 2*7 + 4, 0*7 + 4},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 0, 6*7 + 4, 4*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 4*7 + 4, 4*7 + 6, 2*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,0,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 6*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 2*7 + 4, 0*7 + 4},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 0, 6*7 + 6, 4*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 4*7 + 4, 6*7 + 6, 2*7 + 6},
                                                                                           } ) },

                    // diagonals
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,1,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 6*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 2*7 + 4, 0*7 + 6},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 0, 6*7 + 4, 4*7 + 4},
                                                                                                   { 2*7 + 4, 4*7 + 4, 4*7 + 6, 0*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,0,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 4*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 0, 2*7 + 2, 2*7 + 4, 0*7 + 4},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 2, 6*7 + 6, 4*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 4*7 + 4, 6*7 + 6, 2*7 + 6},
                                                                                           } ) },

                    // three-corners
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,1,0},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 4*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 2*7 + 4, 0*7 + 6},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 2, 6*7 + 4, 4*7 + 4},
                                                                                                   { 2*7 + 4, 4*7 + 4, 4*7 + 6, 0*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {0,1,1,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 4*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 0, 2*7 + 2, 2*7 + 4, 0*7 + 4},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 2, 6*7 + 4, 4*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 4*7 + 4, 4*7 + 6, 2*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,0,1,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 6*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 2*7 + 4, 0*7 + 4},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 0, 6*7 + 4, 4*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 4*7 + 4, 4*7 + 6, 2*7 + 6},
                                                                                                   { 4*7 + 4, 6*7 + 4, 6*7 + 6, 4*7 + 6},
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > ( tableMarkedNodes<4> {1,1,0,1},
                                                                                           {
                                                                                                   { 0*7 + 0, 2*7 + 0, 2*7 + 2, 0*7 + 2},
                                                                                                   { 2*7 + 0, 4*7 + 0, 4*7 + 2, 2*7 + 2},
                                                                                                   { 4*7 + 0, 6*7 + 0, 6*7 + 2, 4*7 + 2},
                                                                                                   { 0*7 + 2, 2*7 + 2, 2*7 + 4, 0*7 + 4},
                                                                                                   { 2*7 + 2, 4*7 + 2, 4*7 + 4, 2*7 + 4},
                                                                                                   { 4*7 + 2, 6*7 + 2, 6*7 + 6, 4*7 + 4},
                                                                                                   { 0*7 + 4, 2*7 + 4, 2*7 + 6, 0*7 + 6},
                                                                                                   { 2*7 + 4, 4*7 + 4, 6*7 + 6, 2*7 + 6},
                                                                                           } ) },
            };

    std::map<tableMarkedNodes<8>, std::array<int, 8>, compare_tableMarkedNodes<8> > Refinement_tables::permutNodesOfRegion_3b_3D =
            {
                    // all
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,0,0,0,0}, {0,1,2,3,4,5,6,7} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,1,1,1,1,1,1,1}, {0,1,2,3,4,5,6,7} ) },

                    // corners
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,0,0,0,0,0,0,0}, {0,1,2,3,4,5,6,7} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,1,0,0,0,0,0,0}, {1,2,3,0,5,6,7,4} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,1,0,0,0,0,0}, {2,3,0,1,6,7,4,5} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,1,0,0,0,0}, {3,0,1,2,7,4,5,6} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,1,0,0,0}, {4,7,6,5,0,3,2,1} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,0,1,0,0}, {5,4,7,6,1,0,3,2} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,0,0,1,0}, {6,5,4,7,2,1,0,3} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,0,0,0,1}, {7,6,5,4,3,2,1,0} ) },

                    // edges
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,1,0,0,0,0,0,0}, {0,1,2,3,4,5,6,7} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,1,1,0,0,0,0,0}, {1,2,3,0,5,6,7,4} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,1,1,0,0,0,0}, {2,3,0,1,6,7,4,5} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,0,0,1,0,0,0,0}, {3,0,1,2,7,4,5,6} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,0,0,0,1,0,0,0}, {4,0,3,7,5,1,2,6} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,1,0,0,0,1,0,0}, {1,5,6,2,0,4,7,3} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,1,0,0,0,1,0}, {2,6,7,3,1,5,4,0} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,1,0,0,0,1}, {3,7,4,0,2,6,5,1} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,1,1,0,0}, {5,4,7,6,1,0,3,2} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,0,1,1,0}, {6,5,4,7,2,1,0,3} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,0,0,1,1}, {7,6,5,4,3,2,1,0} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,1,0,0,1}, {4,7,6,5,0,3,2,1} ) },

                    // faces
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,1,1,1,0,0,0,0}, {0,1,2,3,4,5,6,7} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,0,0,1,1,1,1}, {4,7,6,5,0,3,2,1} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,0,0,1,1,0,0,1}, {0,3,7,4,1,2,6,5} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,1,1,0,0,1,1,0}, {1,5,6,2,0,4,7,3} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {1,1,0,0,1,1,0,0}, {0,4,5,1,3,7,6,2} ) },
                    { std::pair<tableMarkedNodes<8>, std::array<int, 8> > ( tableMarkedNodes<8> {0,0,1,1,0,0,1,1}, {3,2,6,7,0,1,5,4} ) },
            };

    std::map<tableMarkedNodes<8>, int, compare_tableMarkedNodes<8> > Refinement_tables::nbNodesOnRegionFromExternal_3b_3D =
            {
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {0,0,0,0,0,0,0,0},  0) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,0,0,0,0,0,0,0},  6) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,0,0,0,0,0,0}, 16) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,1,1,0,0,0,0}, 36) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,1,1,1,1,1,1}, 48) },
            };

//    std::map<tableMarkedNodes<8>, std::vector<int>, compare_tableMarkedNodes<8> > Refinement_tables::nodesOnRegionFromExternal_3b_3D =
//            {
//
//                    { std::pair<tableMarkedNodes<8>, std::vector<int> >(tableMarkedNodes<8>{0,0,0,0,0,0,0,0},
//                                                                                       {
//
//                                                                                       })},
//                    { std::pair<tableMarkedNodes<8>, std::vector<int> >(tableMarkedNodes<8>{1,0,0,0,0,0,0,0},
//                                                                        {
//                                                                                {
//                                                                                        2*7*7 + 0*7 + 0,
//                                                                                        0*7*7 + 2*7 + 0,
//                                                                                        0*7*7 + 0*7 + 2,
//                                                                                        2*7*7 + 2*7 + 0,
//                                                                                        0*7*7 + 2*7 + 2,
//                                                                                        2*7*7 + 0*7 + 2},
//                                                                        })},
//                    { std::pair<tableMarkedNodes<8>, std::vector<int> >(tableMarkedNodes<8>{1,1,0,0,0,0,0,0},
//                                                                        {
//                                                                                {
//
//                                                                                        2*7*7 + 0*7 + 0,
//                                                                                        0*7*7 + 2*7 + 0,
//                                                                                        0*7*7 + 0*7 + 2,
//                                                                                        2*7*7 + 2*7 + 0,
//                                                                                        0*7*7 + 2*7 + 2,
//                                                                                        2*7*7 + 0*7 + 2},
//                                                                        })},
//            };

    std::map<tableMarkedNodes<8>, int, compare_tableMarkedNodes<8> > Refinement_tables::nbNodesOnRegionToBuild_3b_3D =
            {
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {0,0,0,0,0,0,0,0},  0) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,0,0,0,0,0,0,0},  1) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,0,0,0,0,0,0},  4) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,1,1,0,0,0,0},  8) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,1,1,1,1,1,1},  8) },
            };

    std::array<int, 48> Refinement_tables::nodesOnRegionFromExternal_3b_3D =
            {
                    // edges ordered : bottom, vertical , top
                    2*7*7 + 0*7 + 0,
                    4*7*7 + 0*7 + 0,
                    6*7*7 + 2*7 + 0,
                    6*7*7 + 4*7 + 0,
                    2*7*7 + 6*7 + 0,
                    4*7*7 + 6*7 + 0,
                    0*7*7 + 2*7 + 0,
                    0*7*7 + 4*7 + 0,

                    0*7*7 + 0*7 + 2,
                    0*7*7 + 0*7 + 4,
                    6*7*7 + 0*7 + 2,
                    6*7*7 + 0*7 + 4,
                    6*7*7 + 6*7 + 2,
                    6*7*7 + 6*7 + 4,
                    0*7*7 + 6*7 + 2,
                    0*7*7 + 6*7 + 4,

                    2*7*7 + 0*7 + 6,
                    4*7*7 + 0*7 + 6,
                    6*7*7 + 2*7 + 6,
                    6*7*7 + 4*7 + 6,
                    2*7*7 + 6*7 + 6,
                    4*7*7 + 6*7 + 6,
                    0*7*7 + 2*7 + 6,
                    0*7*7 + 4*7 + 6,

                    // faces ordered : bottom, top, left, right, front, back
                    2*7*7 + 2*7 + 0,
                    4*7*7 + 2*7 + 0,
                    2*7*7 + 4*7 + 0,
                    4*7*7 + 4*7 + 0,

                    2*7*7 + 2*7 + 6,
                    4*7*7 + 2*7 + 6,
                    2*7*7 + 4*7 + 6,
                    4*7*7 + 4*7 + 6,

                    0*7*7 + 2*7 + 2,
                    0*7*7 + 4*7 + 2,
                    0*7*7 + 2*7 + 4,
                    0*7*7 + 4*7 + 4,

                    6*7*7 + 2*7 + 2,
                    6*7*7 + 4*7 + 2,
                    6*7*7 + 2*7 + 4,
                    6*7*7 + 4*7 + 4,

                    2*7*7 + 0*7 + 2,
                    4*7*7 + 0*7 + 2,
                    2*7*7 + 0*7 + 4,
                    4*7*7 + 0*7 + 4,

                    2*7*7 + 6*7 + 2,
                    4*7*7 + 6*7 + 2,
                    2*7*7 + 6*7 + 4,
                    4*7*7 + 6*7 + 4,
            };

    std::map<tableMarkedNodes<8>, std::vector<std::array<int, 3> >, compare_tableMarkedNodes<8> > Refinement_tables::nodesOnRegionToBuild_3b_3D =
            {
                    { std::pair<tableMarkedNodes<8>, std::vector<std::array<int, 3> > > (
                            tableMarkedNodes<8> {0, 0, 0, 0, 0, 0, 0, 0},
                            {

                            })},
                    { std::pair<tableMarkedNodes<8>, std::vector<std::array<int, 3> > > (
                            tableMarkedNodes<8> {1, 1, 1, 1, 1, 1, 1, 1},
                            {
                                    {2,2,2},
                                    {4,2,2},
                                    {2,4,2},
                                    {4,4,2},
                                    {2,2,4},
                                    {4,2,4},
                                    {2,4,4},
                                    {4,4,4},
                            })},
                    { std::pair<tableMarkedNodes<8>, std::vector<std::array<int, 3> > > (
                            tableMarkedNodes<8> {1, 0, 0, 0, 0, 0, 0, 0},
                            {
                                    {2,2,2},
                            })},
                    { std::pair<tableMarkedNodes<8>, std::vector<std::array<int, 3> > > (
                            tableMarkedNodes<8> {1, 1, 0, 0, 0, 0, 0, 0},
                            {
                                    {2,2,2},
                                    {4,2,2},
                                    {2,4,4},
                                    {4,4,4},
                            })},
                    { std::pair<tableMarkedNodes<8>, std::vector<std::array<int, 3> > > (
                            tableMarkedNodes<8> {1, 1, 1, 1, 0, 0, 0, 0},
                            {
                                    {2,2,2},
                                    {4,2,2},
                                    {2,4,2},
                                    {4,4,2},
                                    {2,2,3},
                                    {4,2,3},
                                    {2,4,3},
                                    {4,4,3},
                            })},
            };

    std::map<tableMarkedNodes<8>, int, compare_tableMarkedNodes<8> > Refinement_tables::nbRegionsOnRegionToBuild_3b_3D =
            {
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {0,0,0,0,0,0,0,0},  1) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,0,0,0,0,0,0,0},  4) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,0,0,0,0,0,0}, 11) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,1,1,0,0,0,0}, 22) },
                    { std::pair<tableMarkedNodes<8>, int > ( tableMarkedNodes<8> {1,1,1,1,1,1,1,1}, 27) },
            };

    std::map<tableMarkedNodes<8>, std::vector< std::array<int, 8> > > Refinement_tables::regionsOnRegionToBuild_3b_3D =
            {
                    { std::pair<tableMarkedNodes<8>, std::vector< std::array<int, 8> > > ( tableMarkedNodes<8> {0,0,0,0,0,0,0,0},
                                                                                           {
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           0*7*7 + 6*7 + 6}
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<8>, std::vector< std::array<int, 8> > > ( tableMarkedNodes<8> {1,0,0,0,0,0,0,0},
                                                                                           {
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2}, // 0
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 6*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           2*7*7 + 2*7 + 2}, // 1
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           0*7*7 + 6*7 + 6}, // 2
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           6*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           0*7*7 + 6*7 + 6}, // 3
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<8>, std::vector< std::array<int, 8> > > ( tableMarkedNodes<8> {1,1,0,0,0,0,0,0},
                                                                                           {
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2},  // 0
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           4*7*7 + 0*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2},  // 1
                                                                                                   {
                                                                                                           4*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2},  // 2
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           0*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           0*7*7 + 6*7 + 6},  // 3
                                                                                                   {
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4},  // 4
                                                                                                   {
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           6*7*7 + 2*7 + 0,
                                                                                                           6*7*7 + 6*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           4*7*7 + 4*7 + 4},  // 5
                                                                                                   {
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           6*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 6*7 + 0,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           0*7*7 + 6*7 + 6},  // 6
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           0*7*7 + 6*7 + 6},  // 7
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4},  // 8
                                                                                                   {
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           4*7*7 + 4*7 + 4},  // 9
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           0*7*7 + 6*7 + 6},  // 10
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<8>, std::vector< std::array<int, 8> > > ( tableMarkedNodes<8> {1,1,1,1,0,0,0,0},
                                                                                           {
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2},  // 0
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           4*7*7 + 0*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2},  // 1
                                                                                                   {
                                                                                                           4*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2},  // 2
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           0*7*7 + 4*7 + 0,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           0*7*7 + 4*7 + 2},  // 3
                                                                                                   {
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2},  // 4
                                                                                                   {
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           6*7*7 + 2*7 + 0,
                                                                                                           6*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2},  // 5
                                                                                                   {
                                                                                                           0*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2,
                                                                                                           0*7*7 + 6*7 + 2},  // 6
                                                                                                   {
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 6*7 + 0,
                                                                                                           2*7*7 + 6*7 + 0,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2},  // 7
                                                                                                   {
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           6*7*7 + 4*7 + 0,
                                                                                                           6*7*7 + 6*7 + 0,
                                                                                                           4*7*7 + 6*7 + 0,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 6*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2},  // 8
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           2*7*7 + 2*7 + 3,
                                                                                                           0*7*7 + 2*7 + 4},  // 9
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 2*7 + 3,
                                                                                                           2*7*7 + 2*7 + 3},  // 10
                                                                                                   {
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 2*7 + 3},  // 11
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           0*7*7 + 4*7 + 2,
                                                                                                           0*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 2*7 + 3,
                                                                                                           2*7*7 + 4*7 + 3,
                                                                                                           0*7*7 + 4*7 + 4},  // 12
                                                                                                   {
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 2*7 + 3,
                                                                                                           4*7*7 + 2*7 + 3,
                                                                                                           4*7*7 + 4*7 + 3,
                                                                                                           2*7*7 + 4*7 + 3},  // 13
                                                                                                   {
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 2*7 + 3,
                                                                                                           6*7*7 + 2*7 + 4,
                                                                                                           6*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 4*7 + 3},  // 14
                                                                                                   {
                                                                                                           0*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2,
                                                                                                           0*7*7 + 6*7 + 2,
                                                                                                           0*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 3,
                                                                                                           2*7*7 + 6*7 + 4,
                                                                                                           0*7*7 + 6*7 + 6},  // 15
                                                                                                   {
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2,
                                                                                                           2*7*7 + 4*7 + 3,
                                                                                                           4*7*7 + 4*7 + 3,
                                                                                                           4*7*7 + 6*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4},  // 16
                                                                                                   {
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 6*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2,
                                                                                                           4*7*7 + 4*7 + 3,
                                                                                                           6*7*7 + 4*7 + 4,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           4*7*7 + 6*7 + 4},  // 17
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 2*7 + 3,
                                                                                                           2*7*7 + 4*7 + 3,
                                                                                                           0*7*7 + 4*7 + 4,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4,
                                                                                                           0*7*7 + 6*7 + 6},  // 18
                                                                                                   {
                                                                                                           4*7*7 + 2*7 + 3,
                                                                                                           6*7*7 + 2*7 + 4,
                                                                                                           6*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 4*7 + 3,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           4*7*7 + 6*7 + 4},  // 19
                                                                                                   {
                                                                                                           2*7*7 + 2*7 + 3,
                                                                                                           4*7*7 + 2*7 + 3,
                                                                                                           4*7*7 + 4*7 + 3,
                                                                                                           2*7*7 + 4*7 + 3,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 6*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4},  // 20
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 6*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           0*7*7 + 6*7 + 6},  // 21
                                                                                           } ) },
                    { std::pair<tableMarkedNodes<8>, std::vector< std::array<int, 8> > > ( tableMarkedNodes<8> {1,1,1,1,1,1,1,1},
                                                                                           {
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2},  // 0
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 0,
                                                                                                           4*7*7 + 0*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2},  // 1
                                                                                                   {
                                                                                                           4*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 0*7 + 0,
                                                                                                           6*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2},  // 2
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           0*7*7 + 4*7 + 0,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           0*7*7 + 4*7 + 2},  // 3
                                                                                                   {
                                                                                                           2*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2},  // 4
                                                                                                   {
                                                                                                           4*7*7 + 2*7 + 0,
                                                                                                           6*7*7 + 2*7 + 0,
                                                                                                           6*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2},  // 5
                                                                                                   {
                                                                                                           0*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           2*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 6*7 + 0,
                                                                                                           0*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2,
                                                                                                           0*7*7 + 6*7 + 2},  // 6
                                                                                                   {
                                                                                                           2*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           4*7*7 + 6*7 + 0,
                                                                                                           2*7*7 + 6*7 + 0,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2},  // 7
                                                                                                   {
                                                                                                           4*7*7 + 4*7 + 0,
                                                                                                           6*7*7 + 4*7 + 0,
                                                                                                           6*7*7 + 6*7 + 0,
                                                                                                           4*7*7 + 6*7 + 0,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 6*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2},  // 8
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           0*7*7 + 0*7 + 4,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           2*7*7 + 2*7 + 4,
                                                                                                           0*7*7 + 2*7 + 4},  // 9
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 2*7 + 4},  // 10
                                                                                                   {
                                                                                                           4*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 0*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           6*7*7 + 0*7 + 4,
                                                                                                           6*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 2*7 + 4},  // 11
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           0*7*7 + 4*7 + 2,
                                                                                                           0*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           0*7*7 + 4*7 + 4},  // 12
                                                                                                   {
                                                                                                           2*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4},  // 13
                                                                                                   {
                                                                                                           4*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 2*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 2*7 + 4,
                                                                                                           6*7*7 + 2*7 + 4,
                                                                                                           6*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4},  // 14
                                                                                                   {
                                                                                                           0*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2,
                                                                                                           0*7*7 + 6*7 + 2,
                                                                                                           0*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4,
                                                                                                           0*7*7 + 6*7 + 4},  // 15
                                                                                                   {
                                                                                                           2*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2,
                                                                                                           2*7*7 + 6*7 + 2,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 6*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4},  // 16
                                                                                                   {
                                                                                                           4*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 4*7 + 2,
                                                                                                           6*7*7 + 6*7 + 2,
                                                                                                           4*7*7 + 6*7 + 2,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           6*7*7 + 4*7 + 4,
                                                                                                           6*7*7 + 6*7 + 4,
                                                                                                           4*7*7 + 6*7 + 4},  // 17
                                                                                                   {
                                                                                                           0*7*7 + 0*7 + 4,
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           2*7*7 + 2*7 + 4,
                                                                                                           0*7*7 + 2*7 + 4,
                                                                                                           0*7*7 + 0*7 + 6,
                                                                                                           2*7*7 + 0*7 + 6,
                                                                                                           2*7*7 + 2*7 + 6,
                                                                                                           0*7*7 + 2*7 + 6},  // 18
                                                                                                   {
                                                                                                           2*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           4*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 0*7 + 6,
                                                                                                           4*7*7 + 0*7 + 6,
                                                                                                           4*7*7 + 2*7 + 6,
                                                                                                           2*7*7 + 2*7 + 6},  // 19
                                                                                                   {
                                                                                                           4*7*7 + 0*7 + 4,
                                                                                                           6*7*7 + 0*7 + 4,
                                                                                                           6*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 0*7 + 6,
                                                                                                           6*7*7 + 2*7 + 6,
                                                                                                           4*7*7 + 2*7 + 6},  // 20
                                                                                                   {
                                                                                                           0*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 2*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           0*7*7 + 4*7 + 4,
                                                                                                           0*7*7 + 2*7 + 6,
                                                                                                           2*7*7 + 2*7 + 6,
                                                                                                           2*7*7 + 4*7 + 6,
                                                                                                           0*7*7 + 4*7 + 6},  // 21
                                                                                                   {
                                                                                                           2*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 2*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 2*7 + 6,
                                                                                                           4*7*7 + 2*7 + 6,
                                                                                                           4*7*7 + 4*7 + 6,
                                                                                                           2*7*7 + 4*7 + 6},  // 22
                                                                                                   {
                                                                                                           4*7*7 + 2*7 + 4,
                                                                                                           6*7*7 + 2*7 + 4,
                                                                                                           6*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 2*7 + 6,
                                                                                                           6*7*7 + 2*7 + 6,
                                                                                                           6*7*7 + 4*7 + 6,
                                                                                                           4*7*7 + 4*7 + 6},  // 23
                                                                                                   {
                                                                                                           0*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4,
                                                                                                           0*7*7 + 6*7 + 4,
                                                                                                           0*7*7 + 4*7 + 6,
                                                                                                           2*7*7 + 4*7 + 6,
                                                                                                           2*7*7 + 6*7 + 6,
                                                                                                           0*7*7 + 6*7 + 6},  // 24
                                                                                                   {
                                                                                                           2*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           4*7*7 + 6*7 + 4,
                                                                                                           2*7*7 + 6*7 + 4,
                                                                                                           2*7*7 + 4*7 + 6,
                                                                                                           4*7*7 + 4*7 + 6,
                                                                                                           4*7*7 + 6*7 + 6,
                                                                                                           2*7*7 + 6*7 + 6},  // 25
                                                                                                   {
                                                                                                           4*7*7 + 4*7 + 4,
                                                                                                           6*7*7 + 4*7 + 4,
                                                                                                           6*7*7 + 6*7 + 4,
                                                                                                           4*7*7 + 6*7 + 4,
                                                                                                           4*7*7 + 4*7 + 6,
                                                                                                           6*7*7 + 4*7 + 6,
                                                                                                           6*7*7 + 6*7 + 6,
                                                                                                           4*7*7 + 6*7 + 6},  // 26
                                                                                           } ) },
            };

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
