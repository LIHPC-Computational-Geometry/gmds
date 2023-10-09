/*----------------------------------------------------------------------------*/
#ifndef CLAIRE_PARAMS_H_
#define CLAIRE_PARAMS_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
/*----------------------------------------------------------------------------*/
// gmds File Headers
#include "LIB_GMDS_CLAIRE_export.h"

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
    struct  LIB_GMDS_CLAIRE_API ParamsAero
    {
        enum Dimension
        {
            DIM_2D = 0, /** Maillage de dimension 2 */
            DIM_3D = 1  /** Maillage de dimension 3 */
        };

        enum AlgoPhaseType
        {
            INIT = 0,    		/** Initialisation du maillage              */
            DISTANCE_GEN = 1, 	/** Calcul du champ distance                */
            GRADIENT_GEN = 2, 	/** Calcul du champ gradient                */
            BLOCKING = 3   		/** Conformization via pyramid insertion    */
        };

        /** algo type */
        Dimension dim;
	     /** */
	     int dimension;
        /** phase we start from */
        AlgoPhaseType start_from;
        /** last phase we do */
        AlgoPhaseType stop_at;
        /** flag indicating that we want debug files to be generated */
        bool with_debug_files;
        /** Input mesh file name */
        std::string input_file;
	     /** Output mesh file name */
        std::string output_file;
        /** Output directory */
        std::string output_dir;

	     /** Input 3D surface mesh file name */
	     std::string input_file_3D_surface;		// 3D parameter
	     /** Block Surface 3D */
	     int block_surface_3D;							// 3D parameter

	     /** Epaisseur de couche limite */
	     double delta_cl;
	     /** Angle d'incidence (in degrees) */
	     double angle_attack;

	     /** Nombre mini de blocs */
        int nbrMinBloc;
	     /** Allows or not the insertions on the blocking */
	     bool insertions_allowed=true;
	     /** Limite x amont/aval */
	     double x_lim;
	     /** Limite y amont/aval */
	     double y_lim;
	     /** Limite z amont/aval */
	     double z_lim;
	     /** Nombre de couches lors de l'extrusion */
	     int nbr_couches;
	     /** Allows or not the insertions on the first layer of blocks */
	     bool insertions_allowed_on_first_layer=false;

	     /** Choose the way the vectors field is computed for the extrusion */
	     int vectors_field;
	     /** Choose the x value of the first zone [-inf, x_VectorField_Z1] */
	     double x_VectorField_Z1;
	     /** Choose the x value of the second zone [x_VectorField_Z2, +inf] */
	     double x_VectorField_Z2;

	     /** Block discretization on the wall */
	     double edge_size_wall ;
	     /** Default block discretization in the domain */
	     double edge_size_default ;
	     /** First edge size in the direction ortho to the wall */
	     double edge_size_first_ortho_wall ;
	     /** Number of cells in the boundary layer */
	     int nbrCellsInCL;

	     /** Number of iterations of the Yao Smoother */
	     int nbr_iter_smoothing_yao;
	     /** Damping parameter for the Yao Smoother */
	     double damping_smoothing_yao;

	     /** Limit between inlet and outlet for SU2 writer */
	     double x_lim_SU2_inoutlet;

	     /** If true the result mesh will have a cut on the X axis */
	     bool axisymetry;
    };
/*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* CLAIRE_PARAMS_H_ */
/*----------------------------------------------------------------------------*/
