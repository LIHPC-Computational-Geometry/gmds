//
// Created by rochec on 23/08/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroEllipticSmoothing_2D.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/claire/Utils.h>
#include <gmds/smoothy/EllipticSmoother2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
/*------------------------------------------------------------------------*/

AeroEllipticSmoothing_2D::AeroEllipticSmoothing_2D(Mesh *AMesh, Variable<int>* Alayer, cad::FACManager* Amanager, cad::GeomMeshLinker* Alinker) {
	m_mesh = AMesh;
	m_layer = Alayer;
	m_manager = Amanager;
	m_linker = Alinker;
}


/*------------------------------------------------------------------------*/
AeroEllipticSmoothing_2D::STATUS
AeroEllipticSmoothing_2D::execute()
{
	//==================================================================
	// REORIENT THE FACES
	//==================================================================
	MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();
	for(auto f_id:m_mesh->faces()) {
		Face f=m_mesh->get<Face>(f_id);
		if (f.normal().dot(math::Vector3d({.0, .0, 1.0})) <= 0) {
			std::vector<TCellID> ns = f.getIDs<Node>();
			std::vector<TCellID> ns2(ns.size());
			for (auto i = 0; i < ns.size(); i++)
				ns2[ns.size() - 1 - i] = ns[i];
			f.set<Node>(ns2);
		}
	}


	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH AND THE NODES ON THE
	// FIRST LAYER
	//==================================================================
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator2D op(m_mesh);
	auto mark_node_NAN = m_mesh->newMark<Node>();
	auto mark_node_on_pnt = m_mesh->newMark<Node>();
	auto mark_node_on_crv = m_mesh->newMark<Node>();
	auto mark_edge_on_crv = m_mesh->newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	auto mark_locked_nodes = m_mesh->newMark<Node>();
	for (auto n_id : m_mesh->nodes()) {
		if (m_mesh->isMarked<Node>(n_id, mark_node_on_crv) || m_mesh->isMarked<Node>(n_id, mark_node_on_pnt)) {
			m_mesh->mark<Node>(n_id, mark_locked_nodes);
		}
		if (m_layer->value(n_id) == 1)
		{
			m_mesh->mark<Node>(n_id, mark_locked_nodes);
		}


		// Add the nodes id of the exterior boundary to the slipping nodes list
		if ( (m_mesh->isMarked<Node>(n_id, mark_node_on_crv) || m_mesh->isMarked<Node>(n_id, mark_node_on_pnt))
		    && m_layer->value(n_id) > 1 )
		{
			m_slippingNodesId.push_back(n_id);
		}
	}

	for (int i=0; i<20; i++) {
		//==================================================================
		// PERFORM THE ELLIPTIC SMOOTHING
		//==================================================================
		EllipticSmoother2D smoother2D(m_mesh);
		smoother2D.lock(mark_locked_nodes);
		smoother2D.execute();

		BoundarySlipping();
	}

	return AeroEllipticSmoothing_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroEllipticSmoothing_2D::BoundarySlipping(){

	// Save the old coords
	std::map<TCellID, math::Point> old_cords;
	for (auto n_id:m_slippingNodesId)
	{
		Node n = m_mesh->get<Node>(n_id);
		old_cords[n_id] = n.point();
	}

	// Compute the new positions of the nodes
	for (auto n_id:m_slippingNodesId)
	{
		Node n = m_mesh->get<Node>(n_id);
		std::vector<Node> neighbors_nodes = math::Utils::AdjacentNodes(m_mesh, n) ;
		std::vector<Node> slipping_neighbors_nodes;
		for (auto neighbor_node:neighbors_nodes)
		{
			if (m_layer->value(n_id) == m_layer->value(neighbor_node.id()))
			{
				slipping_neighbors_nodes.push_back(neighbor_node);
			}
		}

		math::Point M = math::Utils::WeightedPointOnBranch(old_cords[slipping_neighbors_nodes[0].id()], n.point(),
		                                                   old_cords[slipping_neighbors_nodes[1].id()], 0.5) ;

		// Projection of M on the curve
		int geom_id = m_linker->getGeomId<Node>(n.id()) ;
		int geom_dim = m_linker->getGeomDim<Node>(n.id()) ;
		if (geom_dim == 2)
		{
			cad::GeomCurve* curve = m_manager->getCurve(geom_id);
			curve->project(M);
			n.setPoint(M);
		}

	}

}
/*------------------------------------------------------------------------*/