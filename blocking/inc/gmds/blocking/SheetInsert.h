/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_SHEETINSERT_H
#define GMDS_BLOCKING_SHEETINSERT_H
/*----------------------------------------------------------------------------*/
#include <map>
#include <string>
/*----------------------------------------------------------------------------*/
//#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
//#include <gmds/math/Point.h>
#include <gmds/blocking/Blocking.h>
#include "LIB_GMDS_BLOCKING_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace blocking{
/*----------------------------------------------------------------------------*/
class LIB_GMDS_BLOCKING_API SheetInsert
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
	SheetInsert();

	/** @brief  Destructor
	 */
	virtual ~SheetInsert();

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS execute();

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS pillow();

	void setBl(gmds::blocking::Blocking* ABl) {bl_ = ABl;};
	gmds::blocking::Blocking* bl() const {return bl_;};
	LCC_3* lcc() {return bl_->lcc();};


 private:
	// blocking structure
	gmds::blocking::Blocking* bl_;

	// shrink set

};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_SHEETINSERT_H
/*----------------------------------------------------------------------------*/


