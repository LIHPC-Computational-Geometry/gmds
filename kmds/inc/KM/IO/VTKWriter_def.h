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