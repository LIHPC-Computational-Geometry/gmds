//
// Created by rochec on 21/03/2022.
//
/*----------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroMeshQuality.h>
#include <gmds/ig/Blocking2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <Eigen/Sparse>
#include <Eigen/Eigen>
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
			std::cout << "Noeud isolé : " << n_id << std::endl;
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

		double ar = AeroMeshQuality::AspectRatioQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double iad = AeroMeshQuality::InternalAngleDeviationQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double eas = AeroMeshQuality::EquiAngleSkewnessQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double cond = AeroMeshQuality::ConditionQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double er = AeroMeshQuality::EdgeRatioQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double jac = AeroMeshQuality::JacobianQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double scajac = AeroMeshQuality::ScaledJacobianQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double shape = AeroMeshQuality::ShapeQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double skew = AeroMeshQuality::SkewQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());
		double stretch = AeroMeshQuality::StretchQUAD(face_nodes[0].point(), face_nodes[1].point(), face_nodes[2].point(), face_nodes[3].point());

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
void Utils::BuildMesh2DFromBlocking2D(Blocking2D* blocking2D, Mesh* m, int mark_block_nodes, int mark_first_layer_nodes, int mark_farfield_nodes){

	m->clear();

	std::map<TCellID, TCellID> map_new_node_ids;
	Variable<int>* var_couche = blocking2D->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche");
	Variable<int>* var_couche_mesh = m->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	// Create all the nodes in the mesh m
	for (auto n_id:blocking2D->nodes())
	{
		Node n_blocking = blocking2D->get<Node>(n_id);
		Node n_mesh = m->newNode(n_blocking.point());
		map_new_node_ids[n_blocking.id()] = n_mesh.id();
		var_couche_mesh->set(n_mesh.id(), var_couche->value(n_id));
		//std::cout << "Node " << n_mesh.id() << ", couche " << var_couche->value(n_id) << ", couche mesh " << var_couche_mesh->value(n_mesh.id()) << std::endl;

		if ( var_couche->value(n_id)==1 || var_couche->value(n_id)==0 )
		{
			m->mark(n_mesh, mark_first_layer_nodes);
		}
		if ( blocking2D->getBlockingDim(n_blocking.id()) == 0 )
		{
			m->mark(n_mesh, mark_block_nodes);
		}

	}

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


	// Mark nodes on the farfield
	// Get max layer id
	int max_layer_id(0);
	for (auto n_id:blocking2D->nodes())
	{
		max_layer_id = std::max(max_layer_id, var_couche->value(n_id));
	}

	for (auto b:blocking2D->allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();

		if ( var_couche->value(b(0,0).id()) == max_layer_id
		    && var_couche->value(b(Nx-1,0).id()) == max_layer_id  )
		{
			for (int i=0; i<Nx; i++)
			{
				m->mark(b(i,0), mark_farfield_nodes);
			}
		}

		if ( var_couche->value(b(0,Ny-1).id()) == max_layer_id
		    && var_couche->value(b(Nx-1,Ny-1).id()) == max_layer_id  )
		{
			for (int i=0; i<Nx; i++)
			{
				m->mark(b(i,Ny-1), mark_farfield_nodes);
			}
		}

		if ( var_couche->value(b(0,0).id()) == max_layer_id
		    && var_couche->value(b(0,Ny-1).id()) == max_layer_id  )
		{
			for (int j=0; j<Ny; j++)
			{
				m->mark(b(0,j), mark_farfield_nodes);
			}
		}

		if ( var_couche->value(b(Nx-1,0).id()) == max_layer_id
		    && var_couche->value(b(Nx-1,Ny-1).id()) == max_layer_id  )
		{
			for (int j=0; j<Ny; j++)
			{
				m->mark(b(Nx-1,j), mark_farfield_nodes);
			}
		}

	}

	gmds::MeshDoctor doc(m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Point Utils::WeightedPointOnBranch(const math::Point A, const math::Point B, const math::Point C, double alpha) {
	math::Point P_Weighted;
	math::Vector3d Vec_AB = B-A ;
	math::Vector3d Vec_BC = C-B ;
	double norme_1 = Vec_AB.norm() ;
	double norme_2 = Vec_BC.norm() ;
	double norme_branche = norme_1 + norme_2 ;
	double norme_cible = alpha*norme_branche ;

	if (norme_cible <= norme_1){
		Vec_AB.normalize();
		P_Weighted = A + norme_cible*Vec_AB ;
	}
	else if (norme_cible > norme_1){
		math::Vector3d Vec_CB = - Vec_BC ;
		Vec_CB.normalize();
		P_Weighted = C + (norme_branche-norme_cible)*Vec_CB ;
	}
	else
	{
		P_Weighted = B ;
	}

	return P_Weighted;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
bool Utils::isInTriangle(const math::Point T1, const math::Point T2, const math::Point T3, const math::Point M)
{
	bool isInFace(false);

	math::Vector3d vij = T2-T1 ;
	math::Vector3d vjk = T3-T2 ;
	math::Vector3d vki = T1-T3 ;
	math::Vector3d viM = M-T1 ;
	math::Vector3d vjM = M-T2 ;
	math::Vector3d vkM = M-T3 ;

	double d1 = ( vij.cross(viM) ).dot( viM.cross(-vki) ) ;
	double d2 = ( -vij.cross(vjM) ).dot( vjM.cross(vjk) ) ;
	double d3 = ( vki.cross(vkM) ).dot( vkM.cross(-vjk) ) ;

	if (d1 >= 0 && d2 >= 0 && d3 >= 0) {
		isInFace = true;
	}
	return isInFace;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double Utils::linearInterpolation2D3Pt(const math::Point P1, const math::Point P2, const math::Point P3, const math::Point M, const double c1, const double c2, const double c3)
{
	Eigen::Matrix3d Mat_A;

	Mat_A(0,0) = P1.X() ;
	Mat_A(0,1) = P1.Y() ;
	Mat_A(0,2) = 1.0 ;
	Mat_A(1,0) = P2.X() ;
	Mat_A(1,1) = P2.Y() ;
	Mat_A(1,2) = 1.0 ;
	Mat_A(2,0) = P3.X() ;
	Mat_A(2,1) = P3.Y() ;
	Mat_A(2,2) = 1.0 ;

	Eigen::Matrix3d Mat_A_Inv = Mat_A.inverse();

	Eigen::Vector3d b;
	b[0] = c1 ;
	b[1] = c2 ;
	b[2] = c3 ;

	Eigen::Vector3d coef = Mat_A_Inv * b;

	return coef[0]*M.X() + coef[1]*M.Y() + coef[2];

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void Utils::CurveBlockEdgesReavel(Blocking2D* blocking2D, Mesh* m){

	m->clear();

	// Create all the faces in the mesh m
	for (auto b:blocking2D->allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();
		for (int i=0; i < Nx-1; i++) {

			Node n0 = m->newNode( b(i,0).point() );
			Node n1 = m->newNode( b(i+1,0).point());
			Face f1 = m->newTriangle(n0, n1, n1);

			Node n2 = m->newNode( b(i,Ny-1).point() );
			Node n3 = m->newNode( b(i+1,Ny-1).point() );
			Face f2 = m->newTriangle(n2, n3, n3);

		}

		for (int j=0; j < Ny-1; j++) {

			Node n0 = m->newNode( b(0,j).point() );
			Node n1 = m->newNode( b(0,j+1).point() );
			Face f1 = m->newTriangle(n0, n1, n1);

			Node n2 = m->newNode( b(Nx-1,j).point() );
			Node n3 = m->newNode(b(Nx-1,j+1).point() );
			Face f2 = m->newTriangle(n2, n3, n3);

		}


	}


}
/*------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/