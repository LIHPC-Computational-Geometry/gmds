
#ifndef GMDS_ELLIPTICMORPH_H
#define GMDS_ELLIPTICMORPH_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MORPHMESH_export.h"
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace morphmesh {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_MORPHMESH_API EllipticMorph
{
	/*
	 * Need a parameter file in argument
	 *
	 * param file format:
	 *
	 * Input file name 				//string
	 * <a_input_filename>.vtk|.mli2
	 *
	 * Output file name 				//string
	 * <a_output_filename>.vtk|.mli2
	 *
	 * To morph							//string
	 * <groupname0>
	 * <groupname1>
	 * ...
	 *
	 * External lock					//string
	 * <groupname2>
	 * <groupname3>
	 * ...
	 *
	 * Internal lock					//string
    * <groupname4>
    * <groupname5>
    * ...
    *
    * Ellipses 						//double double double double
    * X0 coefa0 coefb0 coefc0
	 * X1 coefa1 coefb1 coefc1
	 * ...
	 */

 public:
	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
	 * @param AFilename parameters file name
	 * @param AMesh Empty mesh used to instantiate intern mesh attribute
	 */
	EllipticMorph(const std::string& AFilename, Mesh* AMesh);

	/*------------------------------------------------------------------------*/
	/** @brief Destructor.
	 */
	~EllipticMorph();

	/*------------------------------------------------------------------------*/
	/** @brief Execute the elliptic algorithm using parameters given in argument.
	 */
	void execute();

 private:

	std::vector<TCellID> withExteriorLock();

	std::vector<TCellID> noExteriorLock();

	/*Surement a retirer plus tard*/
	void markLockedCells();

	/** @brief End the elliptic morph algorithm by writing the mesh result in a .vtk or .mli2 file
	 */
	void finalize();

 private:

	/* Mesh mark for locked nodes */
	int m_locked;
	//List of gmds nodes that will be moved in the algorithm
	std::vector<Node> m_morphed;
	//List of gmds nodes that are locked as internal origin
	std::vector<Node> m_lockedIn;
	//List of gmds nodes that are locked as external surface
	std::vector<Node> m_lockedOut;

	std::string m_inputName;
	std::string m_outputName;

	//Lists used to store groups name (extracted from param file) which are used in the algorithm
	//Groups name that will be moved in the algorithm
	std::vector<std::string> m_to_morph;
	//Lists used to store groups name (extracted from param file) which are used in the algorithm
	//Groups name that are locked as external surface
	std::vector<std::string> m_ext_lock;
	//Lists used to store groups name (extracted from param file) which are used in the algorithm
	//Groups name that are locked as internal origin
	std::vector<std::string> m_int_lock;

	/* *Ellipses (extracted from param file) store as lists of doubles as follows {X, coefa, coefb, coefc}, X is the
	 * position of the ellispe on X axis and coefs are
    *                a
    *                |
    *                |
    *     b ------------------- b
    *                |
    *                |
    *                c
	 */
	std::vector<std::vector<double>> m_ellipses;

	std::map<TCellID, math::Vector3d> m_vecs;

	int m_morphRegions;

	int m_lockRegions;
	std::vector<Face> outerSkin;

	std::vector<TCellID> m_not_treated_nodes;

	double m_minEdgeLength;

	/* The mesh to deform */
	Mesh* m_mesh;

};
/*----------------------------------------------------------------------------*/
}// namespace morphmesh
/*----------------------------------------------------------------------------*/
}// namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_ELLIPTICMORPH_H
