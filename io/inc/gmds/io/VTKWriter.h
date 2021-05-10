/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_VTKWRITER_H
#define GMDS_VTKWRITER_H
/*----------------------------------------------------------------------------*/
#include <map>
#include <sstream>
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <gmds/io/IWriter.h>
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

    class VTKWriter : public IWriter {
    public:



        /** @brief Constructor
         *
         * @param AMeshService an implementation of an io service to write data
         * 					   into a mesh
         */
        VTKWriter(IMeshIOService *AMeshService);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor.	*/
        virtual ~VTKWriter();

    protected:

        void initialize(const std::string &AFileName);
        void writeNodes();
        void writeEdges();
        void writeFaces();
        void writeRegions();
        void writeDataNodes();
        void writeDataEdges();
        void writeDataFaces();
        void writeDataRegions();

        void finalize();


        std::map<gmds::TCellID, TInt> m_node_ids_mapping;
        std::map<gmds::TCellID, TInt> m_cell_ids_mapping;

        std::vector<IMeshIOService::EdgeInfo> m_edges_info;
        std::vector<IMeshIOService::CellInfo> m_faces_info;
        std::vector<IMeshIOService::CellInfo> m_regions_info;

        // 0 for edge info, 1 for face info, 2 for region info
        IMeshIOService::DataID m_cell_id[3];
        std::vector<IMeshIOService::DataInt   > m_cell_var_int[3];
        std::vector<IMeshIOService::DataReal  > m_cell_var_real[3];
        std::vector<IMeshIOService::DataVector> m_cell_var_vec[3];

        std::ostringstream m_data_node_stream;

    };
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_VTKWRITER_H
