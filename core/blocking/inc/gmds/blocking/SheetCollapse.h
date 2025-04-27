#ifndef GMDS_SHEETCOLLAPSE_H
#define GMDS_SHEETCOLLAPSE_H
/*----------------------------------------------------------------------------*/
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
//#include <gmds/math/Point.h>
#include <gmds/blocking/Blocking.h>
#include "GMDSBlocking_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
/*----------------------------------------------------------------------------*/
class GMDSBlocking_API SheetCollapse
{
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS,
		NOT_YET_IMPLEMENTED
	} STATUS;

	/** @brief  Default Constructor
	 */
	SheetCollapse();

	/** @brief  Destructor
	 */
	virtual ~SheetCollapse();

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS execute();

	void setDart(Dart_handle ADart) {m_origin = ADart;};
	void setBl(gmds::blocking::Blocking* ABl) {bl_ = ABl;};
	gmds::blocking::Blocking* bl() const {return bl_;};
	LCC_3* lcc() {return bl_->lcc();};

 private:
	/* Reimplementation of erase_marked_darts method of the lcc class, needed to bypass error in Generalized_map.h with topo_unsew(Adart, int) method
	 */
	unsigned int erase(LCC_3::size_type amark);

 private:
	// blocking structure
	gmds::blocking::Blocking* bl_;
	Dart_handle m_origin;

};
}
}

#endif     // GMDS_SHEETCOLLAPSE_H
