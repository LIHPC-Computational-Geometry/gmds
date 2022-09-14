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
#include <gmds/cadfac/FACCurve.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cadfac/FACPoint.h>
#include <gmds/cadfac/FACSurface.h>
#include <gmds/cadfac/FACVolume.h>

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
	STATUS execute(LCC_3::size_type AMark);

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS pillow(LCC_3::size_type AMark);

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	int cleanup_flat_cells(LCC_3::size_type AMark);

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS buildCADfromGrid(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy, int ANz);

	void setBl(gmds::blocking::Blocking* ABl) {bl_ = ABl;};
	gmds::blocking::Blocking* bl() const {return bl_;};
	LCC_3* lcc() {return bl_->lcc();};


 private:

	// blocking structure
	gmds::blocking::Blocking* bl_;

	// shrink set


	// blocking to CAD
	std::map<Dart_handle, gmds::cad::FACPoint*> d2p_;
	std::map<Dart_handle, gmds::cad::FACCurve*> d2c_;
	std::map<Dart_handle, gmds::cad::FACSurface*> d2s_;
	std::map<Dart_handle, gmds::cad::FACVolume*> d2v_;

	gmds::cad::FACManager cad_;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_SHEETINSERT_H
/*----------------------------------------------------------------------------*/


