/*----------------------------------------------------------------------------*/
/** \file    IReader.h
 *  \author  F. LEDOUX
 *  \date    03/17/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IREADER_H_
#define GMDS_IREADER_H_
/*----------------------------------------------------------------------------*/
// GMDS header files
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#include "IMeshIOService.h"
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds {

    /**@class IReader
     *
     * @brief Reader functions to implement
     *
     */
    class IReader {
    public:

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor.	*/
        virtual ~IReader();

        /*------------------------------------------------------------------------*/
        /**@brief By specifying a mesh model, we give the indication about which
         *        cells must be read. For instance N|F|E specify that nodes, edges
         *         and faces must be read
         *
         * @param AModel A mesh model where we are only consider flags N, E,F and R
         */
        void setCellOptions(const MeshModel& AModel);

        /*------------------------------------------------------------------------*/
        /**@brief By specifying a mesh, we give the indication about which data
         *        must be read. For instance N|F|E specify that data relative to
         *        nodes, edges and faces must be read
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
        void read(const std::string &AFileName);

    protected:

        /** @brief Constructor
             *
             * @param AMeshService an implementation of an io service to write data
             * 					   into a mesh
             * @param ACellModel type of cells to import
             * @param ADataModel type of cells we want to import the data from
             */
        IReader(IMeshIOService *AMeshService,
                const MeshModel& ACellModel=MeshModel(DIM3),
                const MeshModel& ADataModel=MeshModel(DIM3));

        /** @brief Once the file opened, this method is called at the beginning of
         * the factory method to give the opportunity to check format error
         *
         * @return
         */
        virtual bool preCheckFormat()=0;

        /**
         * @brief Default behaviour that is the inability to read nodes
         */
        virtual void readNodes();
        /**
         * @brief Default behaviour that is the inability to read edges
         */
        virtual void readEdges();
        /**
         * @brief Default behaviour that is the inability to read faces
         */
        virtual void readFaces();
        /**
         * @brief Default behaviour that is the inability to read regions
         */
        virtual void readRegions();

        /**
         * @brief Default behaviour that is the inability to read nodes data
         */
        virtual void readDataNodes();
        /**
         * @brief Default behaviour that is the inability to read edges data
         */
        virtual void readDataEdges();
        /**
         * @brief Default behaviour that is the inability to read faces data
         */
        virtual void readDataFaces();
        /**
         * @brief Default behaviour that is the inability to read regions data
         */
        virtual void readDataRegions();

        /**@brief Move the stream pointer onto the first occurrency of @AString
         *
         * @param AString the word we look for in the file
         *
         * @return true if we find @AString, false otherwise
         */
        bool moveStreamOntoFirst(const std::string &AString);


    protected:

        /** object giving access to mesh service*/
        IMeshIOService *m_mesh_service;
        /** the file stream */
        std::ifstream* m_stream;

        MeshModel m_cell_model;

        MeshModel m_data_model;
    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_IREADER_H_
/*----------------------------------------------------------------------------*/
