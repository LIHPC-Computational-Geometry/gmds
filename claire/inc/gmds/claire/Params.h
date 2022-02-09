/*----------------------------------------------------------------------------*/
#ifndef CLAIRE_PARAMS_H_
#define CLAIRE_PARAMS_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
/*----------------------------------------------------------------------------*/
// gmds File Headers
/*----------------------------------------------------------------------------*/
/** \brief This file gathers all the parameter structures used in the FHeDo
 *         component.
 */
/*----------------------------------------------------------------------------*/
namespace gmds {

/*------------------------------------------------------------------------*/
/** \struct ParamAero
 *  \brief  structure that gathers Aero global parameters
 */
    struct  ParamsAero
    {
        enum Dimension
        {
            DIM_2D = 0,
            DIM_3D = 1 /** .... */
        };

        enum AlgoPhaseType
        {
            INIT = 0,    /** Frame field generation                  */
            DISTANCE_GEN = 1, /** Frame field smoothing                   */
            GRADIENT_GEN = 2, /** Point field generation                  */
            BLOCKING = 3   	/** Conformization via pyramid insertion    */
        };

        /** algo type */
        Dimension dim;
        /** phase we start from */
        AlgoPhaseType start_from;
        /** last phase we do */
        AlgoPhaseType stop_at;
        /** flag indicating that we want debug files to be generated */
        bool with_debug_files;
        /** Input mesh file name */
        std::string input_file;
        /** Output directory */
        std::string output_dir;
    };
/*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* CLAIRE_PARAMS_H_ */
/*----------------------------------------------------------------------------*/
