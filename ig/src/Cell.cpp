/*----------------------------------------------------------------------------*/
/*
 * Cell.cpp
 *
 *  Created on: 5 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Cell.h>
#include <gmds/ig/Node.h>
#include <gmds/ig/Edge.h>
#include <gmds/ig/Face.h>
#include <gmds/ig/Region.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
Cell::Cell(Mesh* AMesh,const ECellType& AType, const TCellID& AID)
:m_owner(AMesh), m_type(AType), m_id(AID)
{}
/*----------------------------------------------------------------------------*/
TCellID Cell::id() const
{
	return m_id;
}
/*----------------------------------------------------------------------------*/
ECellType  Cell::type() const
{
	return m_type;
}
/*----------------------------------------------------------------------------*/
void
Cell::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	if(GMDS_NODE == type()) {
		throw GMDSException("Cell::computeBoundingBox not available for nodes.");
	}

        minXYZ[0] =  HUGE_VAL;
        minXYZ[1] =  HUGE_VAL;
        minXYZ[2] =  HUGE_VAL;
        maxXYZ[0] = -HUGE_VAL;
        maxXYZ[1] = -HUGE_VAL;
        maxXYZ[2] = -HUGE_VAL;

	std::vector<Node> nodes = get<Node>();

        for(auto & node : nodes) {
                if(node.X() < minXYZ[0]) minXYZ[0] = node.X();
                if(node.Y() < minXYZ[1]) minXYZ[1] = node.Y();
                if(node.Z() < minXYZ[2]) minXYZ[2] = node.Z();
                if(node.X() > maxXYZ[0]) maxXYZ[0] = node.X();
                if(node.Y() > maxXYZ[1]) maxXYZ[1] = node.Y();
                if(node.Z() > maxXYZ[2]) maxXYZ[2] = node.Z();
        }
}
/*------------------------------------------------------------------------*/
void Cell::serializeCellData(std::ostream& AStr) const{
	//the mesh owner pointer is not serialize
	AStr.write((char*)&m_type, sizeof(ECellType));
	AStr.write((char*)&m_id	 , sizeof(TCellID)	);
}
/*------------------------------------------------------------------------*/
void Cell::unserializeCellData(std::istream& AStr) {
	AStr.read((char*)&m_type, sizeof(ECellType)	);
	AStr.read((char*)&m_id	, sizeof(TCellID)	);

//	/* we have to free the indirect connectivity before creating new ones. To
//	 * do that, we remove all the connections (direct and indirect).
//	 */
//	removeAllConnections();
//
//	AStr.read((char*)&adjDirect_  [0],descriptor::sizeDirect  *sizeof(id));
//
//	/* now we create new indirection connections:
//	 * we do not care about which type of connections it is (to nodes, edges,
//	 * faces, or regions. We just allocate new memory for that.
//	 */
//	LIDVectorAllocator<GChunkSize>* allocator =
//					this->m_mesh->getUndefinedSizeConnectivityAllocator();
//	for(int i=0;i<descriptor::sizeIndirect;i++){
//		int nb_ids=0;
//		// value the number of ids stored in the ith connection
//		AStr.read((char*)&nb_ids,sizeof(id));
//		if(nb_ids!=0){
//			id* p = allocator->allocate(nb_ids);
//			/* p[0] stores nb_ids ids */
//
//			AStr.read((char*)&p[1],nb_ids*sizeof(id));
//			adjIndirect_[i]=p;
//		}
//	}
}

/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Node>(std::vector<TCellID>& ACells) const {
	delegateGetNodeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Edge>(std::vector<TCellID>& ACells) const {
	delegateGetEdgeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Face>(std::vector<TCellID>& ACells) const {
	delegateGetFaceIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Region>(std::vector<TCellID>& ACells) const {
	delegateGetRegionIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Node>(std::vector<TCellID>& ACells) const {
        delegateGetAllNodeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Edge>(std::vector<TCellID>& ACells) const {
        delegateGetAllEdgeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Face>(std::vector<TCellID>& ACells) const {
        delegateGetAllFaceIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Region>(std::vector<TCellID>& ACells) const {
        delegateGetAllRegionIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Node>(const std::vector<TCellID>& ACells)  {
	delegateSetNodeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Edge>(const std::vector<TCellID>& ACells)  {
	delegateSetEdgeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Face>(const std::vector<TCellID>& ACells)  {
	delegateSetFaceIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Region>(const std::vector<TCellID>& ACells)  {
	delegateSetRegionIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Node>(const std::vector<Node>& ACells)  {
	delegateSetNodeIDs(convertCellToID<Node>(ACells));
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Edge>(const std::vector<Edge>& ACells)  {
	delegateSetEdgeIDs(convertCellToID<Edge>(ACells));
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Face>(const std::vector<Face>& ACells)  {
	delegateSetFaceIDs(convertCellToID<Face>(ACells));
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Cell::set<Region>(const std::vector<Region>& ACells)  {
	delegateSetRegionIDs(convertCellToID<Region>(ACells));
}

/*----------------------------------------------------------------------------*/
template<> void Cell::add<Node>(TCellID AElt) {
	delegateNodeAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::add<Edge>(TCellID AElt) {
	delegateEdgeAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::add<Face>(TCellID AElt) {
	delegateFaceAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::add<Region>(TCellID AElt) {
	delegateRegionAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Node>(TCellID AElt) {
	delegateNodeRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Edge>(TCellID AElt) {
	delegateEdgeRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Face>(TCellID AElt) {
	delegateFaceRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Region>(TCellID AElt) {
	delegateRegionRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Node>(TCellID AID1, TCellID AID2) {
	delegateNodeReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Edge>(TCellID AID1, TCellID AID2) {
	delegateEdgeReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Face>(TCellID AID1, TCellID AID2) {
	delegateFaceReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Region>(TCellID AID1, TCellID AID2) {
	delegateRegionReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
