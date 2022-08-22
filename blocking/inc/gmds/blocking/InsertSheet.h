/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_INSERTSHEET_H
#define GMDS_BLOCKING_INSERTSHEET_H
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
class LIB_GMDS_BLOCKING_API InsertSheet{

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
	InsertSheet();

	/** @brief  Destructor
	 */
	virtual ~InsertSheet();

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS execute();

	void setBl(gmds::blocking::Blocking* ABl) {bl = ABl;};
	gmds::blocking::Blocking* getBl() const {return bl;};

 private:
	// blocking structure
	gmds::blocking::Blocking* bl;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_INSERTSHEET_H
/*----------------------------------------------------------------------------*/


