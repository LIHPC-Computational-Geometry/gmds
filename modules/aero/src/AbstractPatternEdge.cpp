//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/AbstractPatternEdge.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractPatternEdge::AbstractPatternEdge(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
                                         Mesh *AMeshT, FastLocalize *Afl,
                                         double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  m_mesh(AMesh),
  m_Front(AFront),
  m_e_id(Ae_id),
  m_StructManager(AStructManager),
  m_meshT(AMeshT),
  m_fl(Afl),
  m_dc(dc),
  m_DistanceField(A_DistanceField),
  m_VectorField(A_VectorField)
{

}
/*------------------------------------------------------------------------*/
AbstractPatternEdge::STATUS AbstractPatternEdge::execute()
{
	computeNewHex();
	return AbstractPatternEdge::SUCCESS;
}
/*------------------------------------------------------------------------*/
std::vector<Region>
AbstractPatternEdge::getNewHex()
{
	return m_hex;
}
/*------------------------------------------------------------------------*/