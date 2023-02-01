/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_INPUTMARKEDDARTS_H
#define GMDS_BLOCKING_INPUTMARKEDDARTS_H
/*----------------------------------------------------------------------------*/
//#include <map>
//#include <string>
/*----------------------------------------------------------------------------*/
//#include <CGAL/Generalized_map.h>
//#include <CGAL/Linear_cell_complex_for_generalized_map.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/Blocking.h>
#include "LIB_GMDS_BLOCKING_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace blocking{
/*----------------------------------------------------------------------------*/
class LIB_GMDS_BLOCKING_API InputMarkedDarts
{
 public:

	/** @brief  Default Constructor
	 */
	InputMarkedDarts();

	/** @brief  Destructor
	 */
	virtual ~InputMarkedDarts();
	/** @brief  Dummy call
	 */
	int pillow_mark_first_cells3d(LCC_3* ALcc, LCC_3::size_type AMark, int ANbCells);
	int pillow_mark_cells_3d(LCC_3* ALcc, LCC_3::size_type AMark, std::set<int> ACellset);

	int insertsheet_mark_first_cells3d(LCC_3* ALcc, LCC_3::size_type AMark, int ANbCells);

	int insertsheet_mark_intersect_3d(LCC_3* ALcc, LCC_3::size_type AMark, int ANi, int ANj, int ANk);

	int mark_faceset_regionset_3d(LCC_3* ALcc, LCC_3::size_type AMark, std::set<int> AFaceset, std::set<int> ARegionset) const;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_INPUTMARKEDDARTS_H
/*----------------------------------------------------------------------------*/


