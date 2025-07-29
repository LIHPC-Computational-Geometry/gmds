//
// Created by calderans on 7/2/25.
//

/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLKWRITER_H
#define GMDS_BLKWRITER_H
/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
#include "Blocking3D.h"
#include "LIB_GMDS_AERO_export.h"
#include "gmds/utils/CommonTypes.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
	class LIB_GMDS_AERO_API BLKWriter {

	 public:
	   /** @brief Constructor
         *
         * @param AMeshService an implementation of an io service to write data
         * 					   into a mesh
	    */
	   BLKWriter(Blocking3D* bl);

	   /*------------------------------------------------------------------------*/
	   /** \brief  Destructor.	*/
	   virtual ~BLKWriter();

	   void write(const std::string& AFilename);

	 private:

	   Blocking3D* m_blocking;

	   std::map<int,int> m_node_ids_mapping;
	   std::map<int,int> m_edge_ids_mapping;
	   std::map<int,int> m_face_ids_mapping;
	   std::map<int,int> m_block_ids_mapping;

	};
}
#endif     // GMDS_BLKWRITER_H
