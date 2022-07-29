/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_H
#define GMDS_BLOCKING_H
/*----------------------------------------------------------------------------*/
//#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_BLOCKING_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace blocking{
/*----------------------------------------------------------------------------*/
typedef CGAL::Linear_cell_complex_for_generalized_map<3,3> LCC_3;
typedef LCC_3::Dart_handle                                 Dart_handle;
typedef LCC_3::Point                                       Point;
typedef LCC_3::Vector                                      Vector;
typedef LCC_3::FT                                          FT;
/*----------------------------------------------------------------------------*/
class LIB_GMDS_BLOCKING_API Blocking{

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
	Blocking();

	/** @brief  Destructor
	 */
	virtual ~Blocking();

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS execute();

 private:

	// linear cell complex
	LCC_3 llc;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_H
/*----------------------------------------------------------------------------*/


