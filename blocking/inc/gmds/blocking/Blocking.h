/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_H
#define GMDS_BLOCKING_H
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
#include <gmds/math/Point.h>
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

	/*--------------------------------------------------------------------*/
	/** @brief  Get nbVertices
	 */
	 int nbVertices() const;

	 /*--------------------------------------------------------------------*/
	 /** @brief  Get nbEdges
	  */
	 int nbEdges() const;

	 /*--------------------------------------------------------------------*/
	 /** @brief  Get nbFaces
	  */
	 int nbFaces() const;

	 /*--------------------------------------------------------------------*/
	 /** @brief  Get nbBlocks
	  */
	 int nbBlocks() const;

	 /*--------------------------------------------------------------------*/
	 /** @brief  Create grid of blocks
	  */
	 void createGrid();

	 /*--------------------------------------------------------------------*/
	 /** @brief  Create grid of blocks
	  */
	 void createGrid(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy, int ANz);

	 /*--------------------------------------------------------------------*/
	 /** @brief  Create grid of blocks
	  */
	 void writeMokaFile(std::string AFileName) const;

 private:

	// blocking entities to darts mapping
	std::map<int, Dart_handle> v2d_;
	std::map<int, Dart_handle> e2d_;
	std::map<int, Dart_handle> f2d_;
	std::map<int, Dart_handle> b2d_;

	// linear cell complex
	LCC_3 lcc_;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_H
/*----------------------------------------------------------------------------*/


