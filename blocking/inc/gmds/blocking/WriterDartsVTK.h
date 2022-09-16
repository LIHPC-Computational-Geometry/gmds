/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_WRITERDARTSVTK_H
#define GMDS_BLOCKING_WRITERDARTSVTK_H
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
#include <gmds/blocking/Blocking.h>
#include "LIB_GMDS_BLOCKING_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace blocking{
/*----------------------------------------------------------------------------*/
class LIB_GMDS_BLOCKING_API WriterDartsVTK
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
	WriterDartsVTK();

	/** @brief  Destructor
	 */
	virtual ~WriterDartsVTK();

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	STATUS execute(std::string AFilename);

	/*--------------------------------------------------------------------*/
	/** @brief  Dummy call
	 */
	void setBl(gmds::blocking::Blocking* ABl) {bl_ = ABl;};
	gmds::blocking::Blocking* bl() const {return bl_;};
	LCC_3* lcc() {return bl_->lcc();};


 private:

	// blocking structure
	gmds::blocking::Blocking* bl_;

	double ratio0_ = 0.9;
	double ratio1_ = 0.9;
	double ratio2_ = 0.9;
	double ratio3_ = 0.85;

	double display_alpha1_ = 0.1;
	double display_alpha2_ = 0.4;
	double display_alpha3_ = 0.5;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_WRITERDARTSVTK_H
/*----------------------------------------------------------------------------*/


