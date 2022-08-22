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
	 // TODO is one better than the other ?
//	 int nbVertices() const {return lcc_.number_of_vertex_attributes();};
	 int nbVertices() const {return lcc_.one_dart_per_cell<0>().size();};

	 /*--------------------------------------------------------------------*/
	 /** @brief  Get nbEdges
	  */
	 int nbEdges() const {return lcc_.one_dart_per_cell<1>().size();};

	 /*--------------------------------------------------------------------*/
	 /** @brief  Get nbFaces
	  */
	 int nbFaces() const {return lcc_.one_dart_per_cell<2>().size();};

	 /*--------------------------------------------------------------------*/
	 /** @brief  Get nbBlocks
	  */
	 int nbBlocks() const {return lcc_.one_dart_per_cell<3>().size();};

	 /*--------------------------------------------------------------------*/
	 /** @brief  Create grid of blocks
	  */
	 void createGrid2d();

	 /*--------------------------------------------------------------------*/
	 /** @brief  Create grid of blocks
	  */
	 void createGrid3d();

	 /*--------------------------------------------------------------------*/
	 /** @brief  Create grid of blocks
	  */
	 void createGrid2d(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy);

	 /*--------------------------------------------------------------------*/
	 /** @brief  Create grid of blocks
	  */
	 void createGrid3d(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy, int ANz);

	 /*--------------------------------------------------------------------*/
	 /** @brief  Read the blocking using the vtk file format
	  */
	 void readVTKFile(std::string AFileName);

	 /*--------------------------------------------------------------------*/
	 /** @brief  Write the blocking using the moka file format
	  */
	 void writeMokaFile(std::string AFileName) const;

	 /*--------------------------------------------------------------------*/
	 /** @brief  Write the blocking using the vtk file format
	  */
	 void writeVTKFile(std::string AFileName) const;

	 // TODO read vtk / create grid with holes /
    // TODO fill the maps entities_2_darts
	 // outside cell ?

	 // TODO insert sheet
	 // TODO detect and delete sheet

	 /*--------------------------------------------------------------------*/
	 /** @brief  Insert sheet
	  */
	 void sheetInsert();

	 /*--------------------------------------------------------------------*/
	 /** @brief  Collapse sheet
	  */
	 void sheetCollapse();

	 // TODO later the same on chords



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