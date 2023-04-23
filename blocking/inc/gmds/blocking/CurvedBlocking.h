/*----------------------------------------------------------------------------*/
#ifndef GMDS_CURVED_BLOCKING_H
#define GMDS_CURVED_BLOCKING_H
/*----------------------------------------------------------------------------*/
#include <CGAL/Generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <string>
#include <gmds/utils/Exception.h>
#include <gmds/utils/CommonTypes.h>
#include <LIB_GMDS_BLOCKING_export.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace blocking{
/*----------------------------------------------------------------------------*/
struct ClassificationInfo{
	int dim; // O: point, 1: curve, 2: surface, 3:volume, 4: non classified
	int id; // id of the geometric entity
	ClassificationInfo(const int ADim=4, const int AId=NullID)
	: dim(ADim), id(AId){}
};
/*----------------------------------------------------------------------------*/
/** When we merge two cells,
 */
struct MergeFunctor {
	template<class Cell_attribute>
	void operator()(Cell_attribute& ca1,Cell_attribute& ca2)
	{
		if (ca1.info().dim == ca2.info().dim) {
			// the cells are classifed on the same dim geom entity
			if (ca1.info().id == ca2.info().id) {
				ca1.info().dim = ca2.info().dim;
			}
			else
				throw GMDSException("Classification error!!!");
		}
		else if (ca1.info().dim < ca2.info().dim) {
			// the cells are classifed on the same dim geom entity
			ca1.info().dim = ca1.info().dim;
			ca1.info().id = ca1.info().id;
		}
		else {     // third case: ca1.info().dim>ca2.info().dim
			ca1.info().dim = ca2.info().dim;
			ca1.info().id = ca2.info().id;
		}
	}
};
/*----------------------------------------------------------------------------*/
struct SplitFunctor{
	template<class Cell_attribute>
	void operator()(Cell_attribute& ca1,Cell_attribute& ca2) {
		ca1.info().dim=ca1.info().dim;
		ca1.info().id=ca1.info().id;
		ca2.info()=ca1.info();
	}
};
/*----------------------------------------------------------------------------*/

struct CellData{
	template<class GMap>
	struct Dart_wrapper {
		typedef CGAL::Cell_attribute<GMap, ClassificationInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor> Node_attribute; // A weight
		typedef CGAL::Cell_attribute<GMap, ClassificationInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor> Edge_attribute; // A weight
		typedef CGAL::Cell_attribute<GMap, ClassificationInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor> Face_attribute; // A weight
		typedef CGAL::Cell_attribute<GMap, ClassificationInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor> Region_attribute; // A weight
		typedef std::tuple<Node_attribute, Edge_attribute, Face_attribute, Region_attribute> Attributes;
	};
};
// Definition of my generalized map.
typedef CGAL::Generalized_map<3, CellData> GMap_3;
typedef GMap_3::Dart_handle             DartHandler;

/*----------------------------------------------------------------------------*/
/**@class CurvedBlocking
 * @brief Provide a curved blocking data structure using the 3-G-Map model
 * 		 as described and provided by CGAL.
 * 		 (see https://doc.cgal.org/latest/Generalized_map/index.html)
 */
class LIB_GMDS_BLOCKING_API CurvedBlocking{

 public:
	/** @brief  Default Constructor
	 */
	CurvedBlocking();

	/** @brief  Destructor
	 */
	virtual ~CurvedBlocking();

	/** Create a single hexahedral block in the blocking structure
	 * @return One of the 48 darts of the block
	 */
	DartHandler createHex();

	void sew(DartHandler d1, DartHandler d2){m_gmap.sew<3>(d1,d2);}
	std::string info() const;
	bool isValidTopology() const {return m_gmap.is_valid();}
 private:
	GMap_3 m_gmap;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
} // namespace gmds


/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKING_H
/*----------------------------------------------------------------------------*/