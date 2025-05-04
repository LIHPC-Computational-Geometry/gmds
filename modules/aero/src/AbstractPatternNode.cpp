//
// Created by rochec on 27/04/23.
//

/*------------------------------------------------------------------------*/
#include <gmds/aero/AbstractPatternNode.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractPatternNode::AbstractPatternNode(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                         Mesh *AMeshT, FastLocalize *Afl,
                                         double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  m_mesh(AMesh),
  m_Front(AFront),
  m_n_id(An_id),
  m_StructManager(AStructManager),
  m_meshT(AMeshT),
  m_fl(Afl),
  m_dc(dc),
  m_DistanceField(A_DistanceField),
  m_VectorField(A_VectorField)
{

}
/*------------------------------------------------------------------------*/
AbstractPatternNode::STATUS AbstractPatternNode::execute()
{
	computeNewHex();
	return AbstractPatternNode::SUCCESS;
}
/*------------------------------------------------------------------------*/
std::vector<Region>
AbstractPatternNode::getNewHex()
{
	return m_hex;
}
/*------------------------------------------------------------------------*/