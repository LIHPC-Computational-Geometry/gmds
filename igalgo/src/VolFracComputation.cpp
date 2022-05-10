/*----------------------------------------------------------------------------*/
#include "gmds/igalgo/VolFracComputation.h"

#include "gmds/igalgo/r2d.h"
#include <gmds/math/Quadrilateral.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
#include <string>
#include <vector>
//#include <set>
/*----------------------------------------------------------------------------*/
void gmds::volfraccomputation_2d(gmds::Mesh *AMesh, const gmds::Mesh *AImprintMesh, gmds::Variable<double>* AVolFrac)
{
	// check validity of the inputs
	bool valid_input = true;
	std::string msg("volfraccomputation_2d ");

	gmds::MeshModel model = AMesh->getModel();
	if (!model.has(gmds::F) || !model.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}

	if(AMesh->getNbFaces() != AMesh->getNbQuadrilaterals()) {
		msg += std::string("AMesh should have only quads");
		valid_input = false;
	}
	gmds::MeshModel modelImprint = AImprintMesh->getModel();
	if (!modelImprint.has(gmds::F) || !modelImprint.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}

	if(AImprintMesh->getNbFaces() != AImprintMesh->getNbTriangles()) {
		msg += std::string("AImprintMesh should have only triangles");
		valid_input = false;
   }

	// check mesh orientation
	// TODO
	for(auto f_id: AMesh->faces()) {
		gmds::Face f = AMesh->get<Face>(f_id);
		std::vector<gmds::Node> n = f.get<gmds::Node>();

		gmds::math::Quadrilateral quad(n[0].point(), n[1].point(), n[2].point(), n[3].point());
		double sj = quad.computeScaledJacobian2D();
		std::cout<<"SJ Face "<<f_id<<"("<<n[0].id()<<", "<<n[1].id()<<", "<<n[2].id()<<", "<<n[3].id()<<": "<<sj<<std::endl;
		if(sj < 0) {
			msg += std::string("AMesh has a bad cell.");
			valid_input = false;
			break;
		}
	}

//	for(auto f_id: AImprintMesh->faces()) {
//		gmds::Face f = AImprintMesh->get<Face>(f_id);
//		std::vector<gmds::Node> n = f.get<gmds::Node>();
//
//		gmds::math::Triangle tri(n[0].point(), n[1].point(), n[2].point());
//		double sj = tri.computeScaledJacobian2D();
//		if(sj < 0) {
//			msg += std::string("AImprintMesh has a bad cell.");
//			valid_input = false;
//			break;
//		}
//	}

	if(!valid_input) {
		throw gmds::GMDSException(msg);
	}


	for(auto tri_id: AImprintMesh->faces()) {
		gmds::Face tri = AImprintMesh->get<Face>(tri_id);
		std::vector<gmds::Node> n_tri = tri.get<gmds::Node>();

		r2d_rvec2 vertices[3];
		vertices[0].x = n_tri[0].X();
		vertices[0].y = n_tri[0].Y();
		vertices[1].x = n_tri[2].X();
		vertices[1].y = n_tri[2].Y();
		vertices[2].x = n_tri[1].X();
		vertices[2].y = n_tri[1].Y();
		r2d_int numverts = 3;
		r2d_plane planes[3];
		r2d_poly_faces_from_verts(planes,  vertices, numverts);

		for(auto f_id: AMesh->faces()) {

			gmds::Face f = AMesh->get<Face>(f_id);

			std::vector<gmds::Node> n_quad = f.get<gmds::Node>();
			gmds::TCoord xyz[4][3];
			xyz[0][0] = n_quad[0].X();
			xyz[0][1] = n_quad[0].Y();
			xyz[0][2] = n_quad[0].Z();
			xyz[1][0] = n_quad[1].X();
			xyz[1][1] = n_quad[1].Y();
			xyz[1][2] = n_quad[1].Z();
			xyz[2][0] = n_quad[2].X();
			xyz[2][1] = n_quad[2].Y();
			xyz[2][2] = n_quad[2].Z();
			xyz[3][0] = n_quad[3].X();
			xyz[3][1] = n_quad[3].Y();
			xyz[3][2] = n_quad[3].Z();

			gmds::math::Point pt0(xyz[0][0], xyz[0][1], xyz[0][2]);
			gmds::math::Point pt1(xyz[1][0], xyz[1][1], xyz[1][2]);
			gmds::math::Point pt2(xyz[2][0], xyz[2][1], xyz[2][2]);
			gmds::math::Point pt3(xyz[3][0], xyz[3][1], xyz[3][2]);
			gmds::math::Quadrilateral q(pt0, pt1, pt2, pt3);

			const double volsurf = q.area();

			r2d_poly poly;
			r2d_rvec2 verts[4];

			verts[0].x = xyz[0][0];
			verts[0].y = xyz[0][1];
			verts[1].x = xyz[1][0];
			verts[1].y = xyz[1][1];
			verts[2].x = xyz[2][0];
			verts[2].y = xyz[2][1];
			verts[3].x = xyz[3][0];
			verts[3].y = xyz[3][1];

			r2d_init_poly(&poly, verts, 4);

			r2d_clip(&poly, planes, 3);
			r2d_int POLY_ORDER = 2;
			r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
			r2d_reduce(&poly, om, POLY_ORDER);

			double vf = AVolFrac->value(f_id);
			//std::cout<<"AVANT MAJ \n"<<"f_id "<<f_id<<" volsurf "<<volsurf<<" vf "<<vf<<" om "<<om[0]/volsurf<<std::endl;
			// add but do not forget to divide by cell area
			AVolFrac->set(f_id, vf+om[0]/volsurf);
			double vfN = AVolFrac->value(f_id);
			//std::cout<<"APRES MAJ \n"<<"f_id "<<f_id<<" volsurf "<<volsurf<<" vfN "<<vfN<<" om "<<om[0]/volsurf<<std::endl;
		}
	}
}
/*----------------------------------------------------------------------------*/
//BoundaryExtractor2D::~BoundaryExtractor2D()
//{
//}
///*----------------------------------------------------------------------------*/
//bool BoundaryExtractor2D::isValid() const
//{
//    //check the input first
//    bool valid_input=true;
//    MeshModel model = m_from_mesh->getModel();
//    if (!model.has(F)  || !model.has(N2F) || !model.has(F2N))
//        valid_input= false;
//
//    if(!valid_input)
//        return false;
//
//    //and now, the output
//    model = m_to_mesh->getModel();
//    if (model.has(R))
//        return false;
//    if (!model.has(E) || !model.has(E2N)| !model.has(N2E))
//        return false;
//
//
//    return true;
//}
//
///*----------------------------------------------------------------------------*/
//void BoundaryExtractor2D::setColorOption(Variable<int> *ANodeColor,
//                                         Variable<int> *AEdgeColor,
//                                         Variable<int>* ANodeColorCurv)
//{
//
//    if(!m_to_mesh->hasVariable(GMDS_NODE, ANodeColor->getName())){
//        std::string mess = "The node variable "+ANodeColor->getName();
//        mess +=" is unknown by the target surface mesh";
//        throw GMDSException(mess);
//    }
//    if(!m_to_mesh->hasVariable(GMDS_NODE, ANodeColorCurv->getName())){
//        std::string mess = "The node variable "+ANodeColorCurv->getName();
//        mess +=" is unknown by the target surface mesh";
//        throw GMDSException(mess);
//    }
//    if(!m_to_mesh->hasVariable(GMDS_EDGE, AEdgeColor->getName())){
//        std::string mess = "The edge variable "+AEdgeColor->getName();
//        mess +=" is unknown by the target surface mesh";
//        throw GMDSException(mess);
//    }
//
//    m_color_node_on_pnt = ANodeColor;
//    m_color_edge_on_curv = AEdgeColor;
//    m_color_node_on_curv = ANodeColorCurv;
//    m_with_color = true;
//}
///*----------------------------------------------------------------------------*/
//void BoundaryExtractor2D::setMappings(std::map<TCellID, TCellID> *ANodeMap,
//                                      std::map<TCellID, TCellID> *AEdgeMap,
//                                      std::map<TCellID, TCellID> *ANodeMapInv,
//                                      std::map<TCellID, TCellID> *AEdgeMapInv)
//{
//    m_node_map      = ANodeMap;
//    m_edge_map      = AEdgeMap;
//    m_node_map_inv  = ANodeMapInv;
//    m_edge_map_inv  = AEdgeMapInv;
//    m_with_mapping  = true;
//}
//
///*----------------------------------------------------------------------------*/
//void BoundaryExtractor2D::execute()
//{
//    //the call to the boundary operator is generic in 2D and 3D for many
//    //treatments. So we use her more marks that we could expect first
//    int mark_node_on_curv   = m_from_mesh->newMark<Node>();
//    int mark_node_on_pnt    = m_from_mesh->newMark<Node>();
//    int mark_node_isolated  = m_from_mesh->newMark<Node>();
//
//    int mark_edge_on_curv   = m_from_mesh->newMark<Edge>();
//
//
//    BoundaryOperator2D boundaryOp(m_from_mesh);
//
//    if (!boundaryOp.isValid()) {
//        std::cout << "Invalid model for boundary operations" << std::endl;
//        throw GMDSException("Invalid model for boundary operations");
//    }
//
//    //==================================================================
//    // Mark boundary cells
//    boundaryOp.markCellOnGeometry(mark_edge_on_curv,
//                                  mark_node_on_curv,
//                                  mark_node_on_pnt,
//                                  mark_node_isolated);
//
//    //==================================================================
//    Variable<int>* from_edge_color =
//            m_from_mesh->newVariable<int, GMDS_EDGE>("BoundaryExtractor3D::exec");
//    Variable<int>* from_node_color =
//            m_from_mesh->newVariable<int, GMDS_NODE>("BoundaryExtractor3D::exec");
//
//    if (m_with_color) {
//        // color  edges on curves, nodes on points
//        boundaryOp.colorEdges(mark_edge_on_curv, mark_node_on_pnt,
//                              from_edge_color);
//        boundaryOp.colorNodes(mark_node_on_pnt, from_node_color);
//    }
//
//    //We have now to build the skin mesh and to add fill colors and mappings.
//    std::map<TCellID , TCellID > node_map;
//
//    //first we build nodes, node mapping and color of nodes on pnts
//    for(auto n_id: m_from_mesh->nodes()){
//        if(m_from_mesh->isMarked<Node>(n_id,mark_node_on_curv)||
//           m_from_mesh->isMarked<Node>(n_id,mark_node_on_pnt)){
//            Node from_n = m_from_mesh->get<Node>(n_id);
//            Node to_n   = m_to_mesh->newNode(from_n.point());
//            if(m_with_mapping){
//                (*m_node_map)[n_id] = to_n.id();
//                (*m_node_map_inv)[to_n.id()] = n_id;
//            }
//            else{
//                node_map[n_id] = to_n.id();
//            }
//            if(m_with_color){
//                (*m_color_node_on_pnt)[to_n.id()] = (*from_node_color)[n_id];
//            }
//        }
//    }
//    //then we build edges, edge mapping and color of edges and nodes on curve
//    for(auto e_id: m_from_mesh->edges()){
//        if(m_from_mesh->isMarked<Edge>(e_id,mark_edge_on_curv)){
//            Edge from_e = m_from_mesh->get<Edge>(e_id);
//            std::vector<TCellID> from_ns = from_e.getIDs<Node>();
//            Edge to_e;
//            if(m_with_mapping){
//                to_e = m_to_mesh->newEdge((*m_node_map)[from_ns[0]],
//                                          (*m_node_map)[from_ns[1]]);
//
//                (*m_edge_map)[e_id] = to_e.id();
//                (*m_edge_map_inv)[to_e.id()] = e_id;
//            }
//            else{
//                to_e = m_to_mesh->newEdge(node_map[from_ns[0]],
//                                          node_map[from_ns[1]]);
//            }
//
//            if(m_with_color){
//                (*m_color_edge_on_curv)[to_e.id()] = (*from_edge_color)[e_id];
//                if(m_with_mapping){
//                    for(auto i=0;i<2;i++) {
//                        if (m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_curv) &&
//                            !m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_pnt)) {
//                            (*m_color_node_on_curv)[(*m_node_map)[from_ns[i]]] =
//                                    (*from_edge_color)[e_id];
//                        }
//                    }
//                }
//                else{
//                    for(auto i=0;i<2;i++) {
//                        if (m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_curv) &&
//                            !m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_pnt)) {
//                            (*m_color_node_on_curv)[node_map[from_ns[i]]] =
//                                    (*from_edge_color)[e_id];
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    m_from_mesh->unmarkAll<Node>(mark_node_on_curv);
//    m_from_mesh->unmarkAll<Node>(mark_node_on_pnt);
//    m_from_mesh->unmarkAll<Node>(mark_node_isolated);
//    m_from_mesh->unmarkAll<Edge>(mark_edge_on_curv);
//
//    m_from_mesh->freeMark<Node>(mark_node_on_curv);
//    m_from_mesh->freeMark<Node>(mark_node_on_pnt);
//    m_from_mesh->freeMark<Node>(mark_node_isolated);
//    m_from_mesh->freeMark<Edge>(mark_edge_on_curv);
//
//    m_from_mesh->deleteVariable(GMDS_NODE, from_node_color);
//    m_from_mesh->deleteVariable(GMDS_EDGE, from_edge_color);
