/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IWRITER_H
#define GMDS_IWRITER_H
/*----------------------------------------------------------------------------*/
// GMDS header files
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#include "IMeshIOService.h"
/*----------------------------------------------------------------------------*/
namespace gmds {

    /**@class IWriter
      *
      * @brief writer functions to implement
      */
    class IWriter {

    public:

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor.	*/
        virtual ~IWriter();

        /*------------------------------------------------------------------------*/
        /**@brief By specifying a mesh model, we give the indication about which
         *        cells must be written. For instance N|F|E specify that nodes, edges
         *         and faces must be written
         *
         * @param AModel A mesh model where we are only consider flags N, E,F and R
         */
        void setCellOptions(const MeshModel& AModel);

        /*------------------------------------------------------------------------*/
        /**@brief By specifying a mesh, we give the indication about which data
         *        must be written. For instance N|F|E specify that data relative to
         *        nodes, edges and faces must be written
         *
         * @param AModel A mesh model where we are only consider flags N, E,F and R
         */
        void setDataOptions(const MeshModel& AModel);

        /*------------------------------------------------------------------------*/
        /** @brief  Factory method to read a file and fill in a gmds mesh
         * 			data structure
         *
         * 	@param AFileName the file to read
         */
        void write(const std::string &AFileName);


    protected:

        /** @brief Constructor
             *
             * @param AMeshService an implementation of an io service to read data
             * 					   into a mesh
             * @param ACellModel type of cells to import
             * @param ADataModel type of cells we want to import the data from
             */
        IWriter(IMeshIOService *AMeshService,
                const MeshModel& ACellModel=MeshModel(DIM3|N),
                const MeshModel& ADataModel=MeshModel(DIM3|N));

        virtual void initialize(const std::string &AFileName)=0;

        virtual void writeNodes();

        virtual void writeEdges();

        virtual void writeFaces();

        virtual void writeRegions();

        virtual void writeDataNodes();

        virtual void writeDataEdges();

        virtual void writeDataFaces();

        virtual void writeDataRegions();

        virtual void finalize() = 0;

    protected:

        /** object giving access to mesh service*/
        IMeshIOService *m_mesh_service;
        /** the file stream */
        std::ofstream* m_stream;

        MeshModel m_cell_model;

        MeshModel m_data_model;

    };

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_IWRITER_H
/*----------------------------------------------------------------------------*/
