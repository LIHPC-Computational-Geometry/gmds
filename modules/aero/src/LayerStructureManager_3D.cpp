//
// Created by rochec on 27/04/23.
//

/*------------------------------------------------------------------------*/
#include <gmds/aero/LayerStructureManager_3D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/aero/Front_3D.h>
#include <gmds/aero/AeroException.h>
#include <gmds/io/IGMeshIOService.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

LayerStructureManager_3D::LayerStructureManager_3D(Mesh *AMeshH, Front_3D* AFront, std::map<TCellID, TCellID> &A_map_new_nodes) :
  m_meshH(AMeshH),
  m_Front(AFront)
{
	m_mark_edgesTreated = m_meshH->newMark<Edge>();
	m_mark_facesTreated = m_meshH->newMark<Face>();
	InitFaceStructInfo(A_map_new_nodes);
	InitEdgeStructInfo();
}
/*------------------------------------------------------------------------*/
LayerStructureManager_3D::~LayerStructureManager_3D()
{
	m_meshH->unmarkAll<Edge>(m_mark_edgesTreated);
	m_meshH->freeMark<Edge>(m_mark_edgesTreated);
	m_meshH->unmarkAll<Face>(m_mark_facesTreated);
	m_meshH->freeMark<Face>(m_mark_facesTreated);
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setFaceTreated(TCellID f_id)
{
	m_meshH->mark(m_meshH->get<Face>(f_id), m_mark_facesTreated);
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setEdgeTreated(TCellID e_id)
{
	m_meshH->mark(m_meshH->get<Edge>(e_id), m_mark_edgesTreated);
}
/*------------------------------------------------------------------------*/
bool
LayerStructureManager_3D::isFaceTreated(TCellID f_id)
{
	return m_meshH->isMarked(m_meshH->get<Face>(f_id), m_mark_facesTreated);
}
/*------------------------------------------------------------------------*/
bool
LayerStructureManager_3D::isEdgeTreated(TCellID e_id)
{
	return m_meshH->isMarked(m_meshH->get<Edge>(e_id), m_mark_edgesTreated);
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setFaceNextNode(TCellID f_id, TCellID n_id, TCellID ne_id)
{
	m_FaceInfo[f_id].next_nodes[n_id] = ne_id;
}
/*------------------------------------------------------------------------*/
TCellID
LayerStructureManager_3D::getFaceNextNode(TCellID f_id, TCellID n_id)
{
	return m_FaceInfo[f_id].next_nodes[n_id];
}
/*------------------------------------------------------------------------*/
TCellID
LayerStructureManager_3D::getFaceIdealNextNode(TCellID f_id, TCellID n_id)
{
	return m_FaceInfo[f_id].next_ideal_nodes[n_id];
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setCornerFaceCreated(TCellID e_id, TCellID n_id)
{
	m_EdgeInfo[e_id].CORNER_n_face_created[n_id] = true;
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setEndFaceCreated(TCellID e_id, TCellID n_id)
{
	m_EdgeInfo[e_id].END_n_face_created[n_id] = true;
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setReversalFacesCreated(TCellID e_id, TCellID n_id)
{
	m_EdgeInfo[e_id].REVERSAL_n_faces_created[n_id] = true;
}
/*------------------------------------------------------------------------*/
bool
LayerStructureManager_3D::isCornerFaceCreated(TCellID e_id, TCellID n_id)
{
	return m_EdgeInfo[e_id].CORNER_n_face_created[n_id];
}
/*------------------------------------------------------------------------*/
bool
LayerStructureManager_3D::isEndFaceCreated(TCellID e_id, TCellID n_id)
{
	return m_EdgeInfo[e_id].END_n_face_created[n_id];
}
/*------------------------------------------------------------------------*/
bool
LayerStructureManager_3D::areReversalFacesCreated(TCellID e_id, TCellID n_id)
{
	return m_EdgeInfo[e_id].REVERSAL_n_faces_created[n_id];
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setCornerNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id, TCellID ne_id)
{
	std::pair<TCellID, TCellID> p(f_id, n_id);
	m_EdgeInfo[e_id].CORNER_next_nodes[p] = ne_id ;
}
/*------------------------------------------------------------------------*/
TCellID
LayerStructureManager_3D::getCornerNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id)
{
	std::pair<TCellID, TCellID> p(f_id, n_id);
	return m_EdgeInfo[e_id].CORNER_next_nodes[p] ;
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setNextDiagNode(TCellID e_id, TCellID n_id, TCellID ne_id)
{
	m_EdgeInfo[e_id].diag_next_node[n_id] = ne_id ;
}
/*------------------------------------------------------------------------*/
TCellID
LayerStructureManager_3D::getNextDiagNode(TCellID e_id, TCellID n_id)
{
	return m_EdgeInfo[e_id].diag_next_node[n_id] ;
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setReversalNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id, TCellID ne_id)
{
	std::pair<TCellID, TCellID> p(f_id, n_id);
	m_EdgeInfo[e_id].REVERSAL_adj_nodes[p] = ne_id ;
}
/*------------------------------------------------------------------------*/
TCellID
LayerStructureManager_3D::getReversalNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id)
{
	std::pair<TCellID, TCellID> p(f_id, n_id);
	return m_EdgeInfo[e_id].REVERSAL_adj_nodes[p] ;
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setReversalNextDiagNode(TCellID e_id, TCellID f_id, TCellID n_id, TCellID ne_id)
{
	std::pair<TCellID, TCellID> p(f_id, n_id);
	m_EdgeInfo[e_id].REVERSAL_diag_nodes[p] = ne_id ;
}
/*------------------------------------------------------------------------*/
TCellID
LayerStructureManager_3D::getReversalNextDiagNode(TCellID e_id, TCellID f_id, TCellID n_id)
{
	std::pair<TCellID, TCellID> p(f_id, n_id);
	return m_EdgeInfo[e_id].REVERSAL_diag_nodes[p] ;
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::setReversalNextMediumNode(TCellID e_id, TCellID n_id, TCellID ne_id)
{
	m_EdgeInfo[e_id].REVERSAL_medium_node[n_id] = ne_id ;
}
/*------------------------------------------------------------------------*/
TCellID
LayerStructureManager_3D::getReversalNextMediumNode(TCellID e_id, TCellID n_id)
{
	return m_EdgeInfo[e_id].REVERSAL_medium_node[n_id] ;
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::InitFaceStructInfo(std::map<TCellID, TCellID> map_new_nodes)
{
	for (auto f_id:m_Front->getFaces())
	{
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();

		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[0].id()] = map_new_nodes[f_nodes[0].id()];
		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[1].id()] = map_new_nodes[f_nodes[1].id()];
		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[2].id()] = map_new_nodes[f_nodes[2].id()];
		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[3].id()] = map_new_nodes[f_nodes[3].id()];

		(*this).setFaceNextNode(f_id, f_nodes[0].id(), map_new_nodes[f_nodes[0].id()]);
		(*this).setFaceNextNode(f_id, f_nodes[1].id(), map_new_nodes[f_nodes[1].id()]);
		(*this).setFaceNextNode(f_id, f_nodes[2].id(), map_new_nodes[f_nodes[2].id()]);
		(*this).setFaceNextNode(f_id, f_nodes[3].id(), map_new_nodes[f_nodes[3].id()]);

	}
}
/*------------------------------------------------------------------------*/
void
LayerStructureManager_3D::InitEdgeStructInfo()
{
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	for (auto e_id:m_meshH->edges())
	{
		Edge e = m_meshH->get<Edge>(e_id);
		std::vector<Node> e_nodes = e.get<Node>() ;
		if (var_node_couche_id->value(e_nodes[0].id()) == m_Front->getFrontID()
		    && var_node_couche_id->value(e_nodes[1].id()) == m_Front->getFrontID())
		{
			// So, the edge is on the front.
			if (var_front_edges_classification->value(e_id)==1)
			{
				// Init Corner Edge information
				m_EdgeInfo[e_id].singularity_type = 1;
				m_EdgeInfo[e_id].CORNER_n_face_created[e_nodes[0].id()] = false;
				m_EdgeInfo[e_id].CORNER_n_face_created[e_nodes[1].id()] = false;
			}
			else if (var_front_edges_classification->value(e_id)==2)
			{
				// Init Corner Edge information
				m_EdgeInfo[e_id].singularity_type = 2;
				m_EdgeInfo[e_id].END_n_face_created[e_nodes[0].id()] = false;
				m_EdgeInfo[e_id].END_n_face_created[e_nodes[1].id()] = false;
			}
		}
	}
}