/*----------------------------------------------------------------------------*/
#ifndef FRAME3D_PARAMS_H_
#define FRAME3D_PARAMS_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
/*----------------------------------------------------------------------------*/
// gmds File Headers
#include <gmds/ig/Mesh.h>
#include <gmds/utils/Parameters.h>
#include "LIB_GMDS_FRAME_3D_export.h"
/*----------------------------------------------------------------------------*/
/** \brief This file gathers all the parameter structures used in the FHeDo
 *         component.
 */
/*----------------------------------------------------------------------------*/
namespace gmds {

/*------------------------------------------------------------------------*/
/** \struct ParamsMark
 *  \brief  structure that gathers FHeDo boolean marks used during the
 *          different algorithms
 */
    struct LIB_GMDS_FRAME_3D_API ParamsMark
    {
        /** node mark for nodes classified on surfaces */
        int mark_node_on_surf;
        /** node mark for nodes classified on curves */
        int mark_node_on_curv;
        /** node mark for nodes classified on points*/
        int mark_node_on_pnt;
        /** node mark for isolated nodes (connected to no cells) */
        int mark_node_isolated;
        /** edge mark for edges classified on surfaces */
        int mark_edge_on_surf;
        /** edge mark for edges classified on curves */
        int mark_edge_on_curv;
        /** face mark for faces classified on surfaces */
        int mark_face_on_surf;
        /** node mark for some nodes that must be hard points, i.e. nodes
         *  in the final mesh*/
        int mark_node_hard;
        /** node mark for some nodes where a frame is defined */
        int mark_node_frame;
    };
/*------------------------------------------------------------------------*/
/** \struct ParamsGlobal
 *  \brief  structure that gathers FHeDo global parameters
 */
    struct LIB_GMDS_FRAME_3D_API ParamsGlobal
    {
        enum AlgoChoice
        {
            STABLE_HEX_GENERATION = 0, /** Extract hexes in stable areas       */
            FULL_TET_RECOMBINATION = 1 /** Recombine tet over all the domain   */
        };

        enum AlgoPhaseType
        {
            FF_GEN = 0,    /** Frame field generation                  */
            FF_SMOOTH = 1, /** Frame field smoothing                   */
            PF_GEN = 2,    /** Point field generation                  */
            TET_GEN = 3,   /** Tet remeshing                           */
            HEX_GEN = 4,   /** Hex generation in stable areas          */
            HEX_DOM = 5,   /** Tet 2 hex recombination                 */
            HEX_CONF = 7   /** Conformization via pyramid insertion    */
        };

        /** algo type */
        AlgoChoice algo_choice;
        /** phase we start from */
        AlgoPhaseType start_from;
        /** last phase we do */
        AlgoPhaseType stop_at;
        /** flag indicating that we want debug files to be generated */
        bool with_debug_files;
        /** Input mesh file name */
        std::string mesh_file;
        /** Output directory */
        std::string output_dir;
        /** Quad boundary must be preserved along all the algorithm*/
        bool with_quad_bnd_constraint;
    };
/*------------------------------------------------------------------------*/
/** \struct ParamsFrameField
 *  \brief  structure that gathers parameters of the Frame Field generation
 */
    struct LIB_GMDS_FRAME_3D_API ParamsFrameField
    {
        enum SolverType {
            OPENNL = 0,
            EIGEN = 1
        };
        enum SmoothingType {
            RAY = 0,
            LIU = 1
        };

        /** solver used for building the frame field */
        SolverType solver_type;
        /** use cotangent_weight or not*/
        bool with_cotangent_weights;
        /** use smoothing or not*/
        bool with_smoothing;
        /** type of smoothing algorithm*/
        SmoothingType smoothing_algo;
        /** Number of smoothing iterations*/
        int smoothing_nb_iter;
        /** Convergence criteria smoothing steps*/
        double smoothing_epsilon;
        /** with or without tet mesh adaptation during FF*/
        bool with_mesh_adaptation;
        /** with or without premeshing around volume singularity lines*/
        bool premeshing_sing_lines;
    };

/*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* FRAME3D_PARAMS_H_ */
/*----------------------------------------------------------------------------*/
