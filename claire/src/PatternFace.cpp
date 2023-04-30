//
// Created by rochec on 28/04/23.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternFace.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternFace::PatternFace(Mesh *AMesh, Front_3D *AFront, TCellID Af_id, LayerStructureManager_3D *AStructManager) :
  m_mesh(AMesh),
  m_Front(AFront),
  m_f_id(Af_id),
  m_StructManager(AStructManager)
{

}
/*------------------------------------------------------------------------*/
PatternFace::STATUS PatternFace::execute()
{
	computeNewHex();
	return PatternFace::SUCCESS;
}
/*------------------------------------------------------------------------*/
void
PatternFace::computeNewHex()
{
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	Face f = m_mesh->get<Face>(m_f_id);
	std::vector<Node> nodes = f.get<Node>();

	TCellID n0_id = m_StructManager->getFaceNextNode(m_f_id, nodes[0].id());
	TCellID n1_id = m_StructManager->getFaceNextNode(m_f_id, nodes[1].id());
	TCellID n2_id = m_StructManager->getFaceNextNode(m_f_id, nodes[2].id());
	TCellID n3_id = m_StructManager->getFaceNextNode(m_f_id, nodes[3].id());

	Node n0 = m_mesh->get<Node>(n0_id);
	Node n1 = m_mesh->get<Node>(n1_id);
	Node n2 = m_mesh->get<Node>(n2_id);
	Node n3 = m_mesh->get<Node>(n3_id);
	var_node_couche_id->set(n0_id, m_Front->getFrontID()+1);
	var_node_couche_id->set(n1_id, m_Front->getFrontID()+1);
	var_node_couche_id->set(n2_id, m_Front->getFrontID()+1);
	var_node_couche_id->set(n3_id, m_Front->getFrontID()+1);

	m_hex.push_back( m_mesh->get<Region>( math::Utils::CreateHexaNConnectivities(m_mesh, nodes[0], nodes[1], nodes[2], nodes[3], n0, n1, n2, n3)));

	TCellID f_new_layer_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n0.id(), n1.id(), n2.id(), n3.id());
	var_face_couche_id->set(f_new_layer_id, m_Front->getFrontID()+1);

	m_StructManager->setFaceTreated(m_f_id);
}
/*------------------------------------------------------------------------*/