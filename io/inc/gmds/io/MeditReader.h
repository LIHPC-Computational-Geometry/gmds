/*----------------------------------------------------------------------------*/
/** \file    MeditReader_def.h
 *  \author  F. LEDOUX
 *  \date    09/11/2008
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDIT_READER_H_
#define GMDS_MEDIT_READER_H_
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <gmds/io/IReader.h>
#include "GMDSIo_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {

	class GMDSIo_API MeditReader : public IReader{
	public:


		/** @brief Constructor
         *
         * @param AMeshService an implementation of an io service to write data
         * 					   into a mesh
         * @param AUseMeditLabels use medit labels to define groups of cells
         * 						  in the mesh
         */
		explicit MeditReader(IMeshIOService *AMeshService,
					bool AUseMeditLabels = false);

		/*------------------------------------------------------------------------*/
		/** \brief  Destructor.	*/
		~MeditReader() override;

		/** @brief Set the option flags for using medit labels or not. NOT USED RIGHT NOW
		 *
		 * @param AWithLabels use medit labels to define groups of cells
         * 					  in the mesh

		 */
		void setOptionWithLabel(bool AWithLabels);
		/*------------------------------------------------------------------------*/
		/** \brief  Read the content of the file named fileName and write it in
         *   		m_mesh.
         *
         *   \param fileName 	name of the file to read
         *   \param mask	 	definition of what we read in the file
         *   				 	(for instance, N|F = nodes + faces)
         *   \param useLabels	uses or not the labels of the medit file
         */
		//void read(const std::string &fileName, int mask, bool useLabels = false);

	protected:

		bool preCheckFormat() override;

		void readNodes() override;
		void readEdges() override;
		void readFaces() override;
 		void readRegions() override;

 		void readTriangles();
 		void readQuadrilaterals();

 		void readTetrahedra();
 		void readHexahedra();

		/** mesh dimension */
		int m_mesh_dimension{};

		bool m_with_labels;
/*
		std::map<int, std::string> m_2D_labels;

 std::map<int, std::string> m_3D_labels;
*/
	};
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDIT_READER_H_
/*----------------------------------------------------------------------------*/
