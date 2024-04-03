#ifndef GMDS_MEDIALAXIS_MEDIALAXIS_H
#define GMDS_MEDIALAXIS_MEDIALAXIS_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace medialaxis{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API Medialaxis{

 public:
	/*-------------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;
	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit Medialaxis();
	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~Medialaxis() =default;
	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	Medialaxis::STATUS execute();
};
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // GMDS_MEDIALAXIS_MEDIALAXIS_H