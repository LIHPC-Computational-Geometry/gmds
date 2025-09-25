#ifndef GMDS_QUADFRONT_QUADFRONT_H
#define GMDS_QUADFRONT_QUADFRONT_H
/*----------------------------------------------------------------------------*/
#include "GMDSquadfront_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace quadfront{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class GMDSquadfront_API Quadfront{

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
	explicit Quadfront();
	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~Quadfront() =default;
	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	Quadfront::STATUS execute();
};
/*----------------------------------------------------------------------------*/
}  // end namespace quadfront
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // GMDS_QUADFRONT_QUADFRONT_H