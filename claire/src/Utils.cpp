//
// Created by rochec on 21/03/2022.
//
/*----------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroMeshQuality.h>
#include <gmds/ig/Blocking2D.h>
#include <gmds/ig/MeshDoctor.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {


/*------------------------------------------------------------------------*/
double Utils::distFromNodeIds(Mesh *AMesh, TCellID n0_id, TCellID n1_id){
	Node n0 = AMesh->get<Node>(n0_id);
	Node n1 = AMesh->get<Node>(n1_id);
	math::Point p0 = n0.point();
	math::Point p1 = n1.point();
	math::Vector3d v = p0-p1;
	//std::cout << "n0 : " << n0_id << ", n1 :" << n1_id << std::endl;
	return v.norm();
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
TCellID Utils::CommonEdge(Mesh *AMesh, TCellID n0_id, TCellID n1_id){
	TCellID e_id(NullID);
	Node n0 = AMesh->get<Node>(n0_id);
	std::vector<Edge> adj_edges = n0.get<Edge>();
	for (auto e:adj_edges){
		if (e_id == NullID){
			Node n_opp = e.getOppositeNode(n0);
			if (n_opp.id() == n1_id){
				e_id = e.id();
			}
		}
	}

	return e_id;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void Utils::MeshCleaner(Mesh *AMesh){
	for (auto n_id:AMesh->nodes())
	{
		Node n = AMesh->get<Node>(n_id);
		if (n.get<Face>().empty()) {
			//std::cout << "Noeud isolé : " << n_id << std::endl;
			AMesh->deleteNode(n_id);
		}
	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<Node> Utils::AdjacentNodes(Mesh* m, Node n){

	std::vector<Edge> adjacent_edges = n.get<Edge>() ;
	std::vector<Node> adjacent_nodes;
	for (auto e:adjacent_edges){
		TCellID ne_id = e.getOppositeNodeId(n);
		Node ne = m->get<Node>(ne_id);
		adjacent_nodes.push_back(ne);
	}

	return adjacent_nodes;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void Utils::AnalyseQuadMeshQuality(Mesh* m)
{
	Variable<double>* var_aspect_ratio = m->newVariable<double, GMDS_FACE>("MQ_Aspect_Ratio");
	Variable<double>* var_int_angle_deviation = m->newVariable<double, GMDS_FACE>("MQ_Internal_Angle_Deviation");
	Variable<double>* var_equi_angle_skewness = m->newVariable<double, GMDS_FACE>("MQ_Equi_Angle_Skewness");
	Variable<double>* var_condition = m->newVariable<double, GMDS_FACE>("MQ_Condition");
	Variable<double>* var_edge_ratio = m->newVariable<double, GMDS_FACE>("MQ_EdgeRatio");
	Variable<double>* var_jacobian = m->newVariable<double, GMDS_FACE>("MQ_Jacobian");
	Variable<double>* var_scaled_jacobian = m->newVariable<double, GMDS_FACE>("MQ_ScaledJacobian");
	Variable<double>* var_shape = m->newVariable<double, GMDS_FACE>("MQ_Shape");
	Variable<double>* var_skew = m->newVariable<double, GMDS_FACE>("MQ_Skew");
	Variable<double>* var_stretch = m->newVariable<double, GMDS_FACE>("MQ_Stretch");

	for (auto f_id:m->faces())
	{
		Face f = m->get<Face>(f_id);
		std::vector<Node> face_nodes = f.get<Node>() ;

		double ar = AeroMeshQuality::AspectRatioQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double iad = AeroMeshQuality::InternalAngleDeviationQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double eas = AeroMeshQuality::EquiAngleSkewnessQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double cond = AeroMeshQuality::ConditionQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double er = AeroMeshQuality::EdgeRatioQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double jac = AeroMeshQuality::JacobianQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double scajac = AeroMeshQuality::ScaledJacobianQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double shape = AeroMeshQuality::ShapeQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double skew = AeroMeshQuality::SkewQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());
		double stretch = AeroMeshQuality::StretchQUAD(m, face_nodes[0].id(), face_nodes[1].id(), face_nodes[2].id(), face_nodes[3].id());

		var_aspect_ratio->set(f_id, ar);
		var_int_angle_deviation->set(f_id, iad);
		var_equi_angle_skewness->set(f_id, eas);
		var_condition->set(f_id, cond);
		var_edge_ratio->set(f_id, er);
		var_jacobian->set(f_id, jac);
		var_scaled_jacobian->set(f_id, scajac);
		var_shape->set(f_id, shape);
		var_skew->set(f_id, skew);
		var_stretch->set(f_id, stretch);

	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void Utils::BuildMesh2DFromBlocking2D(Blocking2D* blocking2D, Mesh* m){

	m->clear();

	std::map<TCellID, TCellID> map_new_node_ids;
	int markTreatedNodes = blocking2D->newMark<Node>();

	// Create all the nodes in the mesh m
	for (auto b:blocking2D->allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();
		for (int i=0; i <= Nx-1; i++)
		{
			for (int j=0; j <= Ny-1; j++)
			{
				Node n = b(i,j);
				if (!blocking2D->isMarked(n, markTreatedNodes))
				{
					Node n_new = m->newNode(n.point());
					map_new_node_ids[n.id()] = n_new.id();
					blocking2D->mark(n, markTreatedNodes);
				}
			}
		}
	}

	blocking2D->unmarkAll<Node>(markTreatedNodes);
	blocking2D->freeMark<Node>(markTreatedNodes);

	// Create all the faces in the mesh m
	for (auto b:blocking2D->allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();
		for (int i=0; i < Nx-1; i++)
		{
			for (int j=0; j < Ny-1; j++)
			{
				Node n0_inBlocking = b(i,j);
				Node n1_inBlocking = b(i+1,j);
				Node n2_inBlocking = b(i+1,j+1);
				Node n3_inBlocking = b(i,j+1);

				Node n0_inMesh = m->get<Node>(map_new_node_ids[n0_inBlocking.id()]);
				Node n1_inMesh = m->get<Node>(map_new_node_ids[n1_inBlocking.id()]);
				Node n2_inMesh = m->get<Node>(map_new_node_ids[n2_inBlocking.id()]);
				Node n3_inMesh = m->get<Node>(map_new_node_ids[n3_inBlocking.id()]);

				Face f = m->newQuad(n0_inMesh, n1_inMesh, n2_inMesh, n3_inMesh);
				// Connectivities N->F (x4)
				n0_inMesh.add<Face>(f.id()); // N->F
				n1_inMesh.add<Face>(f.id()); // N->F
				n2_inMesh.add<Face>(f.id()); // N->F
				n3_inMesh.add<Face>(f.id()); // N->F

			}
		}
	}

	gmds::MeshDoctor doc(m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

}
/*------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/