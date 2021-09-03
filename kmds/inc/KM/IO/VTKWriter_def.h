/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/** \file    VTKWriter_def.h
 *  \author  F. LEDOUX
 *  \date    02/24/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_VTKWRITER_DEF_H_
#define KMDS_VTKWRITER_DEF_H_
/*----------------------------------------------------------------------------*/
template<typename TMesh>
class VTKWriter{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *  \param AMesh the mesh we want to write into a file.
     */
	VTKWriter(TMesh& AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~VTKWriter();

    /*------------------------------------------------------------------------*/
    /** \brief  Write the content of mesh_ into the file named AFileName.
     *
     *  \param AFileName name of the output file without its extension
     *  \param AMask	 we check it to find R and/or F to know which type of
     *  				 vtk file to create.
     */
	void write(const std::string& AFileName, const int& AMask);

protected:

	void writeFaces(const std::string& AFileName);
	void writeRegions(const std::string& AFileName);


	void writeNodes(std::ofstream& AOut);
	void writeCellFaces(std::ofstream& AOut);
	void writeCellRegions(std::ofstream& AOut);
	void writeNodesData(std::ofstream& AOut);
	void writeFacesData(std::ofstream& AOut);
	void writeRegionsData(std::ofstream& AOut);

//	void writeVectorSection(
//		std::ofstream& AOut,
//		const std::vector<gmds::math::Vector>& AV);

	  /*------------------------------------------------------------------------*/
//    /** \brief  Method creating a VTK unstructured grid and saving in .vtu file.
//     */
//	void createUnstructuredGrid(const std::string& AFileName);
//
//    /*------------------------------------------------------------------------*/
//    /** \brief  Method creating a VTK poly data and saving in .vtp file.
//     */
//	void createPolyData(const std::string& AFileName);

	/* a mesh */
	TMesh& mesh_;

	/** mesh dimension */
	int mesh_dimension_;

	/** connection between vtk and gmds node ids. */
	std::map<TCellID,TCellID> nodes_connection_;

	/** connection between vtk and gmds face ids. */
        std::map<TCellID,TCellID> faces_connection_;

	/** connection between vtk and gmds region ids. */
        std::map<TCellID,TCellID> regions_connection_;
};
/*----------------------------------------------------------------------------*/
#endif /* KMDS_VTKWRITER_DEF_H_ */
/*----------------------------------------------------------------------------*/