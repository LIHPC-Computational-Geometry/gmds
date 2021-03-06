/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ELLIPTIC_SMOOTHER_2D_H
#define GMDS_ELLIPTIC_SMOOTHER_2D_H
/*----------------------------------------------------------------------------*/
#include <GMDSSmoothy_export.h>
#include <gmds/ig/Mesh.h>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace smoothy{
/*----------------------------------------------------------------------------*/
/** \class EllipticSmoother2D
 *  \brief This class provides ...
 */
/*----------------------------------------------------------------------------*/
class GMDSSmoothy_API EllipticSmoother2D {
 public:

	/** Constructor
	 * @param AMesh 		the mesh to smooth
	 * @param ALockMark  Boolean mark indicating which nodes remain locked
	 */
	EllipticSmoother2D(Mesh* AMesh);

	/** Destructor. Free the internal node mark
	 */
	virtual ~EllipticSmoother2D();

	/** @brief Lock the boundary nodes	 */
	void lockBoundary();

	/** @brief Lock the nodes marked with the mark @p AMark. In practice, the
	 *  mark is internally duplicated
	 * @param AMark node mark
	 * */
	void lock(const int AMark);

		/** @brief Lock the nodes give in @p AV
		 * @param AV a vector of node ids
		 * */
		void lock(const std::vector<TCellID>& AV);

	void execute();

	bool isValid() const;

 private:
	void initLock();
 private:

	/** the mesh we work on */
	Mesh* m_mesh;
	/** Boolean mark indicating which nodes are locked */
	int m_lock;
};
/*----------------------------------------------------------------------------*/
} // end namespace cad
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_ELLIPTIC_SMOOTHER_2D_H
/*----------------------------------------------------------------------------*/
