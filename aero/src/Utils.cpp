//
// Created by rochec on 21/03/2022.
//
/*----------------------------------------------------------------------------*/
#include <gmds/aero/Utils.h>
#include <gmds/aero/AeroMeshQuality.h>
#include <gmds/aero/AdvectedPointRK4_3D.h>
#include <gmds/ig/Blocking2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/math/BezierCurve.h>
#include <gmds/math/BezierHex.h>
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
void Utils::MeshCleaner(Mesh *AMesh)
{
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
void
Utils::BuildMesh3DFromBlocking3D(Blocking3D* Ablocking, Mesh* Am)
{
	std::map<TCellID, TCellID> map_new_node_ids;
	Variable<int>* var_couche = Ablocking->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche");
	Variable<int>* var_couche_block = Ablocking->getOrCreateVariable<int, GMDS_REGION>("GMDS_Layer");
	Variable<int>* var_couche_mesh = Am->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_couche_mesh_hex = Am->getOrCreateVariable<int, GMDS_REGION>("GMDS_Layer");

	// Create all the nodes in the mesh m
	for (auto n_id:Ablocking->nodes())
	{
		Node n_blocking = Ablocking->get<Node>(n_id);
		Node n_mesh = Am->newNode(n_blocking.point());
		map_new_node_ids[n_blocking.id()] = n_mesh.id();
		var_couche_mesh->set(n_mesh.id(), var_couche->value(n_id));
	}

	// Create all the faces in the mesh m
	for (auto b:Ablocking->allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();
		int Nz = b.getNbDiscretizationK();
		for (int i=0; i < Nx-1; i++)
		{
			for (int j=0; j < Ny-1; j++)
			{
				for (int k=0; k<Nz-1; k++)
				{
					Node n0_inBlocking = b(i, j, k);
					Node n1_inBlocking = b(i+1, j, k);
					Node n2_inBlocking = b(i+1, j+1,k);
					Node n3_inBlocking = b(i, j+1,k);
					Node n4_inBlocking = b(i, j, k+1);
					Node n5_inBlocking = b(i+1, j, k+1);
					Node n6_inBlocking = b(i+1, j+1,k+1);
					Node n7_inBlocking = b(i, j+1,k+1);

					Node n0_inMesh = Am->get<Node>(map_new_node_ids[n0_inBlocking.id()]);
					Node n1_inMesh = Am->get<Node>(map_new_node_ids[n1_inBlocking.id()]);
					Node n2_inMesh = Am->get<Node>(map_new_node_ids[n2_inBlocking.id()]);
					Node n3_inMesh = Am->get<Node>(map_new_node_ids[n3_inBlocking.id()]);
					Node n4_inMesh = Am->get<Node>(map_new_node_ids[n4_inBlocking.id()]);
					Node n5_inMesh = Am->get<Node>(map_new_node_ids[n5_inBlocking.id()]);
					Node n6_inMesh = Am->get<Node>(map_new_node_ids[n6_inBlocking.id()]);
					Node n7_inMesh = Am->get<Node>(map_new_node_ids[n7_inBlocking.id()]);

					Region r = Am->newHex(n0_inMesh, n1_inMesh, n2_inMesh, n3_inMesh,
					                      n4_inMesh, n5_inMesh, n6_inMesh, n7_inMesh);
					// Connectivities N->R (x8)
					n0_inMesh.add<Region>(r.id());     // N->R
					n1_inMesh.add<Region>(r.id());     // N->R
					n2_inMesh.add<Region>(r.id());     // N->R
					n3_inMesh.add<Region>(r.id());     // N->R
					n4_inMesh.add<Region>(r.id());     // N->R
					n5_inMesh.add<Region>(r.id());     // N->R
					n6_inMesh.add<Region>(r.id());     // N->R
					n7_inMesh.add<Region>(r.id());     // N->R

					var_couche_mesh_hex->set(r.id(), var_couche_block->value(b.id()));
				}
			}
		}
	}

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
double
Utils::linearInterpolation3D4Pt(const math::Point& P1, const math::Point& P2, const math::Point& P3, const math::Point& P4, const math::Point& M,
                                       double c1, double c2, double c3, double c4)
{
	Eigen::Matrix4d Mat_A;

	// x1 y1 z1 1
	// x2 y2 z2 1
	// x3 y3 z3 1
	// x4 y4 z4 1
	// Remplissage de la matrice A du système Ax=b
	Mat_A(0,0) = P1.X() ;
	Mat_A(0,1) = P1.Y() ;
	Mat_A(0,2) = P1.Z() ;
	Mat_A(0,3) = 1.0 ;
	Mat_A(1,0) = P2.X() ;
	Mat_A(1,1) = P2.Y() ;
	Mat_A(1,2) = P2.Z() ;
	Mat_A(1,3) = 1.0 ;
	Mat_A(2,0) = P3.X() ;
	Mat_A(2,1) = P3.Y() ;
	Mat_A(2,2) = P3.Z() ;
	Mat_A(2,3) = 1.0 ;
	Mat_A(3,0) = P4.X() ;
	Mat_A(3,1) = P4.Y() ;
	Mat_A(3,2) = P4.Z() ;
	Mat_A(3,3) = 1.0 ;

	Eigen::Matrix4d Inv_Mat_A = Mat_A.inverse();

	Eigen::Vector4d b;
	b[0] = c1 ;
	b[1] = c2 ;
	b[2] = c3 ;
	b[3] = c4 ;

	Eigen::Vector4d coef = Inv_Mat_A * b;

	return coef[0]*M.X() + coef[1]*M.Y() + coef[2]*M.Z() + coef[3];

}
/*------------------------------------------------------------------------*/
void Utils::CurveBlockEdgesReveal(Blocking2D* blocking2D, Mesh* m){

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
void Utils::CurveBlockEdgesReveal3D(Blocking3D* Ablocking3D, Mesh* Am, int Asample)
{
	Am->clear();

	// Create all the faces in the mesh m
	for (auto b:Ablocking3D->allBlocks())
	{
		int nb_I = b.getNbDiscretizationI();
		int nb_J = b.getNbDiscretizationJ();
		int nb_K = b.getNbDiscretizationK();

		Array3D<math::Point> ctrl_pts(nb_I,nb_J,nb_K);

		for (auto i=0;i<nb_I;i++)
		{
			for (auto j=0;j<nb_J;j++)
			{
				for (auto k=0;k<nb_K;k++)
				{
					ctrl_pts(i,j,k) = b(i,j,k).point();
				}
			}
		}

		gmds::math::BezierHex bezier_hex(ctrl_pts);

		Array3D<math::Point> pts = bezier_hex.getDiscretization(Asample-1, Asample-1, Asample-1);

		for (auto i=0;i<Asample-1;i++)
		{
			Node n0 = Am->newNode(pts(i,0,0));
			Node n1 = Am->newNode(pts(i+1,0,0));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(i,0,Asample-1));
			n1 = Am->newNode(pts(i+1,0,Asample-1));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(i,Asample-1,Asample-1));
			n1 = Am->newNode(pts(i+1,Asample-1,Asample-1));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(i,Asample-1,0));
			n1 = Am->newNode(pts(i+1,Asample-1,0));
			Am->newTriangle(n0,n1,n1);
		}
		for (auto j=0;j<Asample-1;j++)
		{
			Node n0 = Am->newNode(pts(0,j,0));
			Node n1 = Am->newNode(pts(0,j+1,0));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(0,j,Asample-1));
			n1 = Am->newNode(pts(0,j+1,Asample-1));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(Asample-1,j,Asample-1));
			n1 = Am->newNode(pts(Asample-1,j+1,Asample-1));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(Asample-1,j,0));
			n1 = Am->newNode(pts(Asample-1,j+1,0));
			Am->newTriangle(n0,n1,n1);
		}
		for (auto k=0;k<Asample-1;k++)
		{
			Node n0 = Am->newNode(pts(0,0,k));
			Node n1 = Am->newNode(pts(0,0,k+1));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(0,Asample-1,k));
			n1 = Am->newNode(pts(0,Asample-1,k+1));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(Asample-1,Asample-1,k));
			n1 = Am->newNode(pts(Asample-1,Asample-1,k+1));
			Am->newTriangle(n0,n1,n1);

			n0 = Am->newNode(pts(Asample-1,0,k));
			n1 = Am->newNode(pts(Asample-1,0,k+1));
			Am->newTriangle(n0,n1,n1);
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
void Utils::cutAxiBlocking2D(Mesh *Amesh){
	std::vector<TCellID> faces2remove;
	for(auto f : Amesh->faces()){
		Face face = Amesh->get<Face>(f);

		if(face.center().Y() < 0){
			faces2remove.push_back(f);
		}
	}
	for(auto f : faces2remove){
		Amesh->deleteFace(f);
	}

	std::vector<TCellID> edges2remove;
	for(auto e : Amesh->edges()){
		Edge edge = Amesh->get<Edge>(e);

		if(edge.center().Y() < 0){
			edges2remove.push_back(e);
		}
	}
	for(auto e : edges2remove){
		Amesh->deleteEdge(e);
	}
	std::vector<TCellID> nodes2remove;
	for(auto n : Amesh->nodes()){
		Node node = Amesh->get<Node>(n);

		if(node.Y() < 0){
			nodes2remove.push_back(n);
		}
	}
	for(auto n : nodes2remove){
		Amesh->deleteNode(n);
	}
}
/*----------------------------------------------------------------------------*/
void
Utils::buildEfromFandConnectivies(Mesh *Amesh)
{
	// N->F
	for (auto f_id:Amesh->faces())
	{
		Face f = Amesh->get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();
		for (auto n:f_nodes)
		{
			n.add<Face>(f_id);
		}
	}

	// Build E
	for (auto f_id:Amesh->faces())
	{
		Face f = Amesh->get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();
		TCellID e_0_id = math::Utils::CommonEdge(Amesh, f_nodes[0].id(), f_nodes[1].id());
		TCellID e_1_id = math::Utils::CommonEdge(Amesh, f_nodes[1].id(), f_nodes[2].id());
		TCellID e_2_id = math::Utils::CommonEdge(Amesh, f_nodes[2].id(), f_nodes[3].id());
		TCellID e_3_id = math::Utils::CommonEdge(Amesh, f_nodes[3].id(), f_nodes[0].id());

		if (e_0_id == NullID)
		{
			Edge e = Amesh->newEdge(f_nodes[0], f_nodes[1]);
			e_0_id = e.id();
			f_nodes[0].add<Edge>(e_0_id);
			f_nodes[1].add<Edge>(e_0_id);
		}
		Edge e_0 = Amesh->get<Edge>(e_0_id);
		e_0.add<Face>(f_id);
		f.add<Edge>(e_0);

		if (e_1_id == NullID)
		{
			Edge e = Amesh->newEdge(f_nodes[1], f_nodes[2]);
			e_1_id = e.id();
			f_nodes[1].add<Edge>(e_1_id);
			f_nodes[2].add<Edge>(e_1_id);
		}
		Edge e_1 = Amesh->get<Edge>(e_1_id);
		e_1.add<Face>(f_id);
		f.add<Edge>(e_1);

		if (e_2_id == NullID)
		{
			Edge e = Amesh->newEdge(f_nodes[2], f_nodes[3]);
			e_2_id = e.id();
			f_nodes[2].add<Edge>(e_2_id);
			f_nodes[3].add<Edge>(e_2_id);
		}
		Edge e_2 = Amesh->get<Edge>(e_2_id);
		e_2.add<Face>(f_id);
		f.add<Edge>(e_2);

		if (e_3_id == NullID)
		{
			Edge e = Amesh->newEdge(f_nodes[3], f_nodes[0]);
			e_3_id = e.id();
			f_nodes[3].add<Edge>(e_3_id);
			f_nodes[0].add<Edge>(e_3_id);
		}
		Edge e_3 = Amesh->get<Edge>(e_3_id);
		e_3.add<Face>(f_id);
		f.add<Edge>(e_3);

	}
}
/*----------------------------------------------------------------------------*/
void
Utils::orientRegion(Mesh *m, Region r)
{
	std::vector<Node> r_nodes = r.get<Node>();
	math::Point r_center = r.center();

	if (r.type()==GMDS_HEX)	// Hexa
	{
		// Check the orientation of one quad
		math::Vector3d v1 = r_nodes[1].point() - r_nodes[0].point() ;
		math::Vector3d v2 = r_nodes[2].point() - r_nodes[0].point() ;

		math::Point p = (r_nodes[0].point() + r_nodes[1].point() + r_nodes[2].point()) ;
		p.setX(p.X()/3.0);
		p.setY(p.Y()/3.0);
		p.setZ(p.Z()/3.0);
		double orient = (v1.cross(v2)).dot(p-r_center);
		if (orient >= 0)	// Normales entrantes je crois ?
		{
			//std::cout << "Hexa reoriented 1" << std::endl;
			std::vector<Node> r_nodes_new;
			r_nodes_new.push_back(r_nodes[0]);
			r_nodes_new.push_back(r_nodes[3]);
			r_nodes_new.push_back(r_nodes[2]);
			r_nodes_new.push_back(r_nodes[1]);
			r_nodes_new.push_back(r_nodes[4]);
			r_nodes_new.push_back(r_nodes[7]);
			r_nodes_new.push_back(r_nodes[6]);
			r_nodes_new.push_back(r_nodes[5]);
			r.set(r_nodes_new);
		}

		// Other faces
		r_nodes = r.get<Node>();
		p = (r_nodes[0].point() + r_nodes[1].point() + r_nodes[5].point()) ;
		p.setX(p.X()/3.0);
		p.setY(p.Y()/3.0);
		p.setZ(p.Z()/3.0);
		v1 = r_nodes[1].point() - r_nodes[0].point() ;
		v2 = r_nodes[5].point() - r_nodes[0].point() ;
		orient = (v1.cross(v2)).dot(p-r_center);
		if (orient >= 0)	// Normales entrantes je crois
		{
			//std::cout << "Hexa reoriented 2" << std::endl;
			std::vector<Node> r_nodes_new;
			r_nodes_new.push_back(r_nodes[4]);
			r_nodes_new.push_back(r_nodes[5]);
			r_nodes_new.push_back(r_nodes[6]);
			r_nodes_new.push_back(r_nodes[7]);
			r_nodes_new.push_back(r_nodes[0]);
			r_nodes_new.push_back(r_nodes[1]);
			r_nodes_new.push_back(r_nodes[2]);
			r_nodes_new.push_back(r_nodes[3]);
			r.set(r_nodes_new);
		}

	}
	else if (r.type()==GMDS_PYRAMID)
	{
		// Check the orientation of the base
		math::Vector3d v1 = r_nodes[1].point() - r_nodes[0].point() ;
		math::Vector3d v2 = r_nodes[2].point() - r_nodes[0].point() ;

		math::Point p = (r_nodes[0].point() + r_nodes[1].point() + r_nodes[2].point()) ;
		p.setX(p.X()/3.0);
		p.setY(p.Y()/3.0);
		p.setZ(p.Z()/3.0);
		double orient = (v1.cross(v2)).dot(p-r_center);
		if (orient <= 0)	// Normales entrantes
		{
			//std::cout << "Pyramid reoriented" << std::endl;
			std::vector<Node> r_nodes_new;
			r_nodes_new.push_back(r_nodes[0]);
			r_nodes_new.push_back(r_nodes[3]);
			r_nodes_new.push_back(r_nodes[2]);
			r_nodes_new.push_back(r_nodes[1]);
			r_nodes_new.push_back(r_nodes[4]);
			r.set(r_nodes_new);
		}
	}
	/*
	else
	{
		std::cout << "WARNING in orientRegion: Region type not supported." << std::endl;
		std::cout << "Type: " << r.type() << std::endl;
	}
	 */
}
/*----------------------------------------------------------------------------*/
double
Utils::BinomialCoefficient(int An, int Ak)
{
	double factorial_n(1);
	double factorial_k(1);
	double factorial_n_k(1);

	for (int i=2; i <= An; i++)
	{
		factorial_n = factorial_n*i;
		if (i <= Ak)
		{
			factorial_k = factorial_k*i;
		}
		if (i <= An-Ak)
		{
			factorial_n_k = factorial_n_k*i;
		}
	}

	return factorial_n/(factorial_k*factorial_n_k);

}
/*------------------------------------------------------------------------*/
double
Utils::BernsteinPolynomial(int An, int Ai, double Au)
{
	double C_n_i = BinomialCoefficient(An, Ai);
	return C_n_i* pow(Au, Ai)* pow(1-Au, An-Ai);
}
/*----------------------------------------------------------------------------*/
double
Utils::DerivativeBernsteinPolynomial(int An, int Ai, double Au)
{
	if (Ai == 0)
	{
		return -An*BernsteinPolynomial(An-1, 0, Au);
	}
	else if (Ai==An)
	{
		return An*BernsteinPolynomial(An-1, An-1, Au);
	}
	else
	{
		return An*(BernsteinPolynomial(An-1, Ai-1, Au)-BernsteinPolynomial(An-1, Ai, Au));
	}

}
/*----------------------------------------------------------------------------*/
math::Vector3d
Utils::DerivativeBezierCurve(std::vector<math::Point> &Pts, double Au)
{
	math::Point deriv;
	for (int i=0;i<Pts.size();i++)
	{
		deriv = deriv + DerivativeBernsteinPolynomial(Pts.size()-1, i, Au)*Pts[i] ;
	}
	math::Vector3d dev;
	dev.setXYZ(deriv.X(), deriv.Y(), deriv.Z());
	return dev;
}
/*----------------------------------------------------------------------------*/
void
Utils::UpdateLinker3D(cad::GeomMeshLinker* linker_1, const Node& n_1, cad::GeomMeshLinker* linker_2, const Node& n_2)
{
	int geom_dim = linker_1->getGeomDim<Node>(n_1.id());
	if(geom_dim == 1)
	{
		linker_2->linkNodeToPoint(n_2.id(), linker_1->getGeomId<Node>(n_1.id()));
	}
	else if(geom_dim==2)
	{
		linker_2->linkNodeToCurve(n_2.id(), linker_1->getGeomId<Node>(n_1.id()));
	}
	else if(geom_dim==3)
	{
		linker_2->linkNodeToSurface(n_2.id(), linker_1->getGeomId<Node>(n_1.id()));
	}
}
/*----------------------------------------------------------------------------*/
void
Utils::resizeMesh(Mesh* Amesh, double Ascale)
{
	for (auto n_id:Amesh->nodes())
	{
		Node n = Amesh->get<Node>(n_id);
		n.setX(n.X()*Ascale);
		n.setY(n.Y()*Ascale);
		n.setZ(n.Z()*Ascale);
	}
}
/*----------------------------------------------------------------------------*/
double
Utils::lengthBezierCurve(BezierCurve* Abc)
{
	int Asample_size(101);
	std::vector<math::Point> sample = Abc->getDiscretization(Asample_size-1);
	double length(0);
	for (int i=0;i<Asample_size-1;i++)
	{
		length += (sample[i+1]-sample[i]).norm();
	}
	return length;
}
/*----------------------------------------------------------------------------*/
double
Utils::maxErrorBtwBezierCurveandGeomCurve(BezierCurve* Abc, cad::GeomCurve* Acurve, int Asample)
{
	double max_error(0.0);
	for (int i=0;i<Asample;i++)
	{
		double param = double(i)/double(Asample-1.0);
		math::Point p_bc = (*Abc)(param) ;
		math::Point p_curve = p_bc;
		Acurve->project(p_curve);
		if (abs( (p_bc-p_curve).norm() ) > max_error)
		{
			max_error = abs( (p_bc-p_curve).norm() );
		}
	}
	return max_error;
}
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/