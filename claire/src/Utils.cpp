//
// Created by rochec on 21/03/2022.
//
/*----------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroMeshQuality.h>
#include <gmds/claire/AdvectedPointRK4_3D.h>
#include <gmds/ig/Blocking2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <Eigen/Sparse>
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
TCellID Utils::CommonEdge(Mesh *AMesh, TCellID n0_id, TCellID n1_id){
	TCellID e_id(NullID);
	Node n0 = AMesh->get<Node>(n0_id);
	std::vector<Edge> adj_edges = n0.get<Edge>();
	for (auto const& e:adj_edges){
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
TCellID Utils::CommonFace(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id){
	TCellID f_id(NullID);
	Node n0 = AMesh->get<Node>(n0_id);
	std::vector<Face> adj_faces = n0.get<Face>();
	for (auto const& f:adj_faces)
	{
		if (f_id==NullID)
		{
			std::vector<Node> f_nodes = f.get<Node>();
			if ( f_nodes.size() == 4
			    && (f_nodes[0].id() == n0_id || f_nodes[0].id() == n1_id || f_nodes[0].id() == n2_id || f_nodes[0].id() == n3_id )
			    && (f_nodes[1].id() == n0_id || f_nodes[1].id() == n1_id || f_nodes[1].id() == n2_id || f_nodes[1].id() == n3_id )
			    && (f_nodes[2].id() == n0_id || f_nodes[2].id() == n1_id || f_nodes[2].id() == n2_id || f_nodes[2].id() == n3_id )
			    && (f_nodes[3].id() == n0_id || f_nodes[3].id() == n1_id || f_nodes[3].id() == n2_id || f_nodes[3].id() == n3_id ))
			{
				f_id = f.id();
			}
		}
	}
	return f_id;

}
/*------------------------------------------------------------------------*/
TCellID Utils::CommonFace3Nodes(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id)
{
	TCellID f_id(NullID);
	Node n0 = AMesh->get<Node>(n0_id);
	std::vector<Face> adj_faces = n0.get<Face>();
	if (adj_faces.size() < 3 || adj_faces.size() > 4)
	{
		std::cout << "Utils::CommonFace3Nodes: This method can't handle other than quad or tri faces." << std::endl;
	}
	for (auto const& f:adj_faces)
	{
		if (f_id==NullID)
		{
			std::vector<Node> f_nodes = f.get<Node>();
			bool condition_quad = ( f_nodes.size() == 4
			    && (f_nodes[0].id() == n0_id || f_nodes[1].id() == n0_id || f_nodes[2].id() == n0_id || f_nodes[3].id() == n0_id )
			    && (f_nodes[0].id() == n1_id || f_nodes[1].id() == n1_id || f_nodes[2].id() == n1_id || f_nodes[3].id() == n1_id )
			    && (f_nodes[0].id() == n2_id || f_nodes[1].id() == n2_id || f_nodes[2].id() == n2_id || f_nodes[3].id() == n2_id ) );
			bool condition_triangle  = ( f_nodes.size() == 3
			    && (f_nodes[0].id() == n0_id || f_nodes[1].id() == n0_id || f_nodes[2].id() == n0_id )
			    && (f_nodes[0].id() == n1_id || f_nodes[1].id() == n1_id || f_nodes[2].id() == n1_id )
			    && (f_nodes[0].id() == n2_id || f_nodes[1].id() == n2_id || f_nodes[2].id() == n2_id ));

			if (condition_quad || condition_triangle)
			{
				f_id = f.id();
			}
		}
	}
	return f_id;
}
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
std::vector<Node> Utils::AdjacentNodes(Mesh* m, const Node& n){

	std::vector<Edge> adjacent_edges = n.get<Edge>() ;
	std::vector<Node> adjacent_nodes;
	for (auto const& e:adjacent_edges){
		TCellID ne_id = e.getOppositeNodeId(n);
		Node ne = m->get<Node>(ne_id);
		adjacent_nodes.push_back(ne);
	}

	return adjacent_nodes;
}
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
void Utils::BuildMesh2DFromBlocking2D(Blocking2D* blocking2D, Mesh* m, TInt mark_block_nodes, TInt mark_first_layer_nodes, TInt mark_farfield_nodes){

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
math::Point Utils::WeightedPointOnBranch(const math::Point& A, const math::Point& B, const math::Point& C, double alpha) {
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
bool Utils::isInTriangle(const math::Point& T1, const math::Point& T2, const math::Point& T3, const math::Point& M)
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
bool Utils::sameSide(const math::Point& T1, const math::Point& T2, const math::Point& T3, const math::Point& M1, const math::Point& M2)
{
	math::Vector3d V1 = T2-T1 ;
	math::Vector3d V2 = T3-T1 ;
	math::Vector3d V3 = M1-T1 ;
	math::Vector3d V4 = M2-T1 ;

	math::Vector3d Normal = V1.cross(V2);
	double dotN4 = Normal.dot(V3);
	double dotM = Normal.dot(V4);

	return (std::signbit(dotN4) == std::signbit(dotM) );
}
/*------------------------------------------------------------------------*/
bool Utils::isInTetra(const math::Point& T1, const math::Point& T2, const math::Point& T3, const math::Point& T4, const math::Point& M)
{
	bool isInTet = false;
	/*
	= (math::Utils::sameSide(T1, T2, T3, T4, M) &&
	                math::Utils::sameSide(T2, T3, T4, T1, M) &&
	                math::Utils::sameSide(T3, T4, T1, T2, M) &&
	                math::Utils::sameSide(T4, T1, T2, T3, M));
	                */

	if (!isInTet)
	{
		TCoord c1, c2, c3, c4;
		Point::computeBarycentric(T1, T2, T3, T4, M, c1, c2, c3, c4);
		//std::cout << c1 << ", " << c2 << ", " << c3 << ", " << c4 << std::endl;
		if (c1 >= 0 && c2 >= 0 && c3 >= 0 && c4 >= 0)
		{
			isInTet = true;
		}
	}

	if (!isInTet)
	{
		isInTet = (math::Utils::sameSide(T1, T2, T3, T4, M) &&
	                math::Utils::sameSide(T2, T3, T4, T1, M) &&
	                math::Utils::sameSide(T3, T4, T1, T2, M) &&
	                math::Utils::sameSide(T4, T1, T2, T3, M));
	}

	/*
	if (!isInTet)
	{
		isInTet = isInTriangle(T1, T2, T3, M);
	}
	if (!isInTet)
	{
		isInTet = isInTriangle(T2, T3, T4, M);
	}
	if (!isInTet)
	{
		isInTet = isInTriangle(T1, T2, T4, M);
	}
	if (!isInTet)
	{
		isInTet = isInTriangle(T1, T3, T4, M);
	}
	 */

	// Test méthode avec les volumes
	/*
	if (!isInTet)
	{
		double V = (1.0/6.0)*abs( T1.X()*(T2.Y()*T3.Z()-T3.Y()*T2.Z())
		                             + T2.X()*( T3.Y()*T1.Z()-T1.Y()*T3.Z())
		                             + T3.X()*(T1.Y()*T2.Z()-T2.Y()*T1.Z())
		                             + T4.X()*(T1.Y()*(T2.Z()-T3.Z())+T2.Y()*(T3.Z()-T1.Z())+T3.Y()*(T1.Z()-T2.Z())) );
		double V1 = (1.0/6.0)*abs( M.X()*(T2.Y()*T3.Z()-T3.Y()*T2.Z())
		                              + T2.X()*( T3.Y()*M.Z()-M.Y()*T3.Z())
		                              + T3.X()*(M.Y()*T2.Z()-T2.Y()*M.Z()) );
		double V2 = (1.0/6.0)*abs( T1.X()*(M.Y()*T3.Z()-T3.Y()*M.Z())
		                              + M.X()*( T3.Y()*T1.Z()-T1.Y()*T3.Z())
		                              + T3.X()*(T1.Y()*M.Z()-M.Y()*T1.Z()) ) ;
		double V3 = (1.0/6.0)*abs( T1.X()*(T2.Y()*M.Z()-M.Y()*T2.Z())
		                              + T2.X()*( M.Y()*T1.Z()-T1.Y()*M.Z())
		                              + M.X()*(T1.Y()*T2.Z()-T2.Y()*T1.Z()) );
		double V4 = (1.0/6.0)*abs(T1.Y()*(M.Y()*(T2.Z()-T3.Z())+T2.Y()*(T3.Z()-M.Z())+T3.Y()*(M.Z()-T2.Z())) );

		if ( abs(V-V1-V2-V3-V4) <= pow(10,-6) )
		{
			isInTet = true;
		}
	}
	 */

	return isInTet;
}
/*------------------------------------------------------------------------*/
double Utils::linearInterpolation2D3Pt(const math::Point& P1, const math::Point& P2, const math::Point& P3, const math::Point& M, double c1, double c2, double c3)
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
TCellID Utils::GetOrCreateEdgeAndConnectivitiesN2E(Mesh *AMesh, TCellID n0_id, TCellID n1_id)
{
	Node n0 = AMesh->get<Node>(n0_id);
	Node n1 = AMesh->get<Node>(n1_id);
	TCellID e_id = math::Utils::CommonEdge(AMesh, n0_id, n1_id);
	if (e_id==NullID)		// Then, the edge doesn't exist yet
	{
		Edge e = AMesh->newEdge(n0, n1);	// E->N (x2)
		n0.add<Edge>(e);	// N->E
		n1.add<Edge>(e);	// N->E
		e_id = e.id();
	}
	return e_id;
}
/*------------------------------------------------------------------------*/
TCellID Utils::GetOrCreateQuadAndConnectivities(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id)
{
	Node n0 = AMesh->get<Node>(n0_id);
	Node n1 = AMesh->get<Node>(n1_id);
	Node n2 = AMesh->get<Node>(n2_id);
	Node n3 = AMesh->get<Node>(n3_id);

	TCellID f_id = math::Utils::CommonFace(AMesh, n0_id, n1_id, n2_id, n3_id);
	if (f_id == NullID)		// Then, the face doesn't exist yet
	{
		Face f = AMesh->newQuad(n0, n1, n2, n3);	// F->N (x4)
		n0.add<Face>(f);	// N->F
		n1.add<Face>(f);	// N->F
		n2.add<Face>(f);	// N->F
		n3.add<Face>(f);	// N->F
		f_id = f.id();

		TCellID e0_id = Utils::GetOrCreateEdgeAndConnectivitiesN2E(AMesh, n0.id(), n1.id());	// E<->N
		TCellID e1_id = Utils::GetOrCreateEdgeAndConnectivitiesN2E(AMesh, n1.id(), n2.id());	// E<->N
		TCellID e2_id = Utils::GetOrCreateEdgeAndConnectivitiesN2E(AMesh, n2.id(), n3.id());	// E<->N
		TCellID e3_id = Utils::GetOrCreateEdgeAndConnectivitiesN2E(AMesh, n3.id(), n0.id());	// E<->N

		Edge e0 = AMesh->get<Edge>(e0_id);
		Edge e1 = AMesh->get<Edge>(e1_id);
		Edge e2 = AMesh->get<Edge>(e2_id);
		Edge e3 = AMesh->get<Edge>(e3_id);

		f.add<Edge>(e0);	// F->E
		f.add<Edge>(e1);	// F->E
		f.add<Edge>(e2);	// F->E
		f.add<Edge>(e3);	// F->E

		e0.add<Face>(f);	// E->F
		e1.add<Face>(f);	// E->F
		e2.add<Face>(f);	// E->F
		e3.add<Face>(f);	// E->F

	}

	return f_id;
}
/*------------------------------------------------------------------------*/
TCellID Utils::CreateHexaNConnectivities(Mesh *Amesh, Node n0, Node n1, Node n2, Node n3, Node n4, Node n5, Node n6, Node n7)
{
	// Create the hex associated to the face
	Region r = Amesh->newHex(n0, n1, n2, n3, n4, n5, n6, n7);	// R->N (x8)
	n0.add<Region>(r);			// N->R (1/8)
	n1.add<Region>(r);			// N->R
	n2.add<Region>(r);			// N->R
	n3.add<Region>(r);			// N->R
	n4.add<Region>(r);			// N->R
	n5.add<Region>(r);			// N->R
	n6.add<Region>(r);			// N->R
	n7.add<Region>(r);			// N->R (8/8)

	// Create the 6 faces and connectivities N<->F (1x F->N and 4x N->F per face)
	TCellID f0_id = math::Utils::GetOrCreateQuadAndConnectivities(Amesh, n0.id(), n1.id(), n2.id(), n3.id());
	TCellID f1_id = math::Utils::GetOrCreateQuadAndConnectivities(Amesh, n4.id(), n5.id(), n6.id(), n7.id());
	TCellID f2_id = math::Utils::GetOrCreateQuadAndConnectivities(Amesh, n0.id(), n1.id(), n5.id(), n4.id());
	TCellID f3_id = math::Utils::GetOrCreateQuadAndConnectivities(Amesh, n1.id(), n2.id(), n6.id(), n5.id());
	TCellID f4_id = math::Utils::GetOrCreateQuadAndConnectivities(Amesh, n2.id(), n3.id(), n7.id(), n6.id());
	TCellID f5_id = math::Utils::GetOrCreateQuadAndConnectivities(Amesh, n3.id(), n0.id(), n4.id(), n7.id());

	Face f0 = Amesh->get<Face>(f0_id);
	Face f1 = Amesh->get<Face>(f1_id);
	Face f2 = Amesh->get<Face>(f2_id);
	Face f3 = Amesh->get<Face>(f3_id);
	Face f4 = Amesh->get<Face>(f4_id);
	Face f5 = Amesh->get<Face>(f5_id);


	// Get the 12 edges
	// Remark: these edges has been created when the 6 faces were created before
	TCellID e0_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n0.id(), n1.id());
	TCellID e1_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n1.id(), n2.id());
	TCellID e2_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n2.id(), n3.id());
	TCellID e3_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n3.id(), n0.id());
	TCellID e4_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n4.id(), n5.id());
	TCellID e5_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n5.id(), n6.id());
	TCellID e6_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n6.id(), n7.id());
	TCellID e7_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n7.id(), n4.id());
	TCellID e8_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n0.id(), n4.id());
	TCellID e9_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n1.id(), n5.id());
	TCellID e10_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n2.id(), n6.id());
	TCellID e11_id = math::Utils::GetOrCreateEdgeAndConnectivitiesN2E(Amesh, n3.id(), n7.id());

	Edge e0 = Amesh->get<Edge>(e0_id);
	Edge e1 = Amesh->get<Edge>(e1_id);
	Edge e2 = Amesh->get<Edge>(e2_id);
	Edge e3 = Amesh->get<Edge>(e3_id);
	Edge e4 = Amesh->get<Edge>(e4_id);
	Edge e5 = Amesh->get<Edge>(e5_id);
	Edge e6 = Amesh->get<Edge>(e6_id);
	Edge e7 = Amesh->get<Edge>(e7_id);
	Edge e8 = Amesh->get<Edge>(e8_id);
	Edge e9 = Amesh->get<Edge>(e9_id);
	Edge e10 = Amesh->get<Edge>(e10_id);
	Edge e11 = Amesh->get<Edge>(e11_id);


	//----------------------------//
	// Create the connectivities 	//
	// 		R<->F x6					//
	//----------------------------//

	f0.add<Region>(r);	// F->R
	r.add<Face>(f0);		// R->F

	f1.add<Region>(r);	// F->R
	r.add<Face>(f1);		// R->F

	f2.add<Region>(r);	// F->R
	r.add<Face>(f2);		// R->F

	f3.add<Region>(r);	// F->R
	r.add<Face>(f3);		// R->F

	f4.add<Region>(r);	// F->R
	r.add<Face>(f4);		// R->F

	f5.add<Region>(r);	// F->R
	r.add<Face>(f5);		// R->F


	//----------------------------//
	// Create the connectivities 	//
	// 		R<->E x12				//
	//----------------------------//

	e0.add<Region>(r);	// E->R
	r.add<Edge>(e0);		// R->E

	e1.add<Region>(r);	// E->R
	r.add<Edge>(e1);		// R->E

	e2.add<Region>(r);	// E->R
	r.add<Edge>(e2);		// R->E

	e3.add<Region>(r);	// E->R
	r.add<Edge>(e3);		// R->E

	e4.add<Region>(r);	// E->R
	r.add<Edge>(e4);		// R->E

	e5.add<Region>(r);	// E->R
	r.add<Edge>(e5);		// R->E

	e6.add<Region>(r);	// E->R
	r.add<Edge>(e6);		// R->E

	e7.add<Region>(r);	// E->R
	r.add<Edge>(e7);		// R->E

	e8.add<Region>(r);	// E->R
	r.add<Edge>(e8);		// R->E

	e9.add<Region>(r);	// E->R
	r.add<Edge>(e9);		// R->E

	e10.add<Region>(r);	// E->R
	r.add<Edge>(e10);		// R->E

	e11.add<Region>(r);	// E->R
	r.add<Edge>(e11);		// R->E

	return r.id();

}
/*------------------------------------------------------------------------*/
bool Utils::isThisHexMeshValid(Mesh *Amesh)
{
	bool isValid(true);

	for (auto r_id:Amesh->regions())
	{
		Region r = Amesh->get<Region>(r_id);
		if (r.get<Node>().size() != 8)	// So the region r is not an HEXA
		{
			isValid = false;
		}
	}

	for (auto f_id:Amesh->faces())
	{
		Face f = Amesh->get<Face>(f_id);
		if (f.get<Node>().size() != 4)	// So the face f is not a QUAD
		{
			isValid = false;
		}
		if (f.area() < 0)
		{
			isValid = false;
		}
	}


	return isValid;
}
/*------------------------------------------------------------------------*/
math::Point Utils::AdvectedPointRK4_UniqVector_3D(Mesh *Amesh, FastLocalize *fl, math::Point M, double dist_cible, Variable<double>* A_distance, math::Vector3d A_vector)
{
	Variable<math::Vector3d>* var_vectors = Amesh->newVariable<math::Vector3d, GMDS_NODE>("Temp_vectors");
	for (auto n_id:Amesh->nodes())
	{
		var_vectors->set(n_id, A_vector);
	}

	AdvectedPointRK4_3D advpoint(Amesh, fl, M, dist_cible, A_distance, var_vectors);
	advpoint.execute();

	Amesh->deleteVariable(GMDS_NODE, var_vectors);

	return advpoint.getPend();
}
/*----------------------------------------------------------------------------*/
std::vector<Face> Utils::getFacesAdjToEdgeInHexa(Mesh *Amesh, TCellID e_id, TCellID r_id)
{
	std::vector<Face> adj_faces;	// Should be sized 2 at the end of the method

	Region r = Amesh->get<Region>(r_id);
	std::vector<Face> r_faces = r.get<Face>();	// Should be sized 6

	if (r_faces.size() != 6)
	{
		std::cout << "ATTENTION getFacesAdjToEdgeInHexa: la région n'est pas un hexa." << std::endl;
		std::cout << "-> Nbr faces: " << r_faces.size() << std::endl;
	}

	for (auto const& f:r_faces)
	{
		std::vector<Edge> f_edges = f.get<Edge>();
		if (f_edges[0].id() == e_id || f_edges[1].id() == e_id
		    || f_edges[2].id() == e_id || f_edges[3].id() == e_id)
		{
			adj_faces.push_back(f);
		}
	}

	if (adj_faces.size() !=2)
	{
		std::cout << "ATTENTION getFacesAdjToEdgeInHexa: les 2 faces adjacentes n'ont pas été trouvées correctement." << std::endl;
	}

	return adj_faces;
}
/*----------------------------------------------------------------------------*/
Edge Utils::oppositeEdgeInFace(Mesh *Amesh, TCellID e_id, TCellID f_id)
{
	Edge e_opp;
	Face f = Amesh->get<Face>(f_id);
	std::vector<Edge> f_edges = f.get<Edge>();
	if (f_edges.size() !=4)
	{
		std::cout << "ATTENTION oppositeEdgeInFace: la face n'est pas un quad." << std::endl;
	}

	Edge e = Amesh->get<Edge>(e_id);
	std::vector<Node> e_nodes = e.get<Node>();

	for (auto const& e_loc:f.get<Edge>())
	{
		std::vector<Node> e_loc_nodes = e_loc.get<Node>();
		if (e_loc_nodes[0].id() != e_nodes[0].id()
		    && e_loc_nodes[0].id() != e_nodes[1].id()
		    && e_loc_nodes[1].id() != e_nodes[0].id()
		    && e_loc_nodes[1].id() != e_nodes[1].id())
		{
			e_opp = e_loc;
		}
	}

	return e_opp;
}
/*----------------------------------------------------------------------------*/
double Utils::minEdgeLenght(Mesh *Amesh){

	double minLenght(std::numeric_limits<double>::max());
	for (auto edge_id:Amesh->edges())
	{
		Edge edge = Amesh->get<Edge>(edge_id);
		if(edge.length() < minLenght)
		{
			minLenght = edge.length() ;
		}
	}
	return minLenght;
}
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/