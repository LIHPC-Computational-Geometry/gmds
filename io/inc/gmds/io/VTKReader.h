/*----------------------------------------------------------------------------*/
/** \file    VTKReader.h
 *  \author  F. LEDOUX
 *  \date    12/17/2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_VTKREADER_H_
#define GMDS_VTKREADER_H_
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <gmds/io/IReader.h>
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

	class VTKReader : public IReader {
	public:



        /** @brief Constructor
         *
         * @param AMeshService an implementation of an io service to write data
         * 					   into a mesh
         */
        VTKReader(IMeshIOService *AMeshService);

        /*------------------------------------------------------------------------*/
		/** \brief  Destructor.	*/
		virtual ~VTKReader();

	protected:

		virtual bool preCheckFormat();

        void readNodes();
        void readEdges();
        void readFaces();
        void readRegions();

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


	private:
	    std::vector<ECellType> m_cell_types;

	    /**Nb max of nodes per cell, which can be either a face or region */
        static const int m_nb_max_node_per_cell;
	};
}
/*----------------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------------*/
