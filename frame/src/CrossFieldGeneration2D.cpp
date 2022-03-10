/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryOperator2D.h>
//#include <GMDS/Algo/DistanceFieldBuilder2D.h>
/*---------------------------------------------------------------------------*/
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <sstream>
/*---------------------------------------------------------------------------*/
#include "gmds/frame/CrossFieldGeneration2D.h"
//#include "StabilityBallCross2D.h"
#include "gmds/frame/LaplaceCross2D.h"
//#include "gmds/frame/FunctionalApproachCross2D.h"
//#include "gmds/frame/LevelSetCross2D.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
namespace {
class FrameLogger
{
 public:
	FrameLogger(bool logEnable) : m_silent(!logEnable) {}
	enum class LogEnd {EndLine, Flush, None};
	void log(const std::string &message, LogEnd logend = LogEnd::EndLine)
	{
		if (!m_silent) {
			switch (logend) {
			case LogEnd::EndLine: std::cout << message << std::endl; break;
			case LogEnd::Flush:   std::cout << message << std::flush; break;
			case LogEnd::None:    std::cout << message; break;
			default: break; }
		}
	}

 private:
	bool m_silent;
};
}     // namespace
/*---------------------------------------------------------------------------*/
CrossFieldGeneration2D::CrossFieldGeneration2D(Mesh* AMesh)
  : m_mesh(AMesh)
{
  if (m_mesh->hasVariable(GMDS_NODE, "cross")){
	m_cross_field_2D = m_mesh->getVariable<math::Cross2D,GMDS_NODE>("cross");
	m_cross_field_2D->clear();
  }
  else {
    m_cross_field_2D = m_mesh->newVariable<math::Cross2D,GMDS_NODE>("cross");
  }

  m_debug_output = "";
}
/*---------------------------------------------------------------------------*/
 void CrossFieldGeneration2D::setDebugEnable(bool ADebugEnable)
{
  m_debugEnable = ADebugEnable;
}/*---------------------------------------------------------------------------*/
 void CrossFieldGeneration2D::setLogEnable(bool ALogEnable)
{
  m_logEnable = ALogEnable;
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::setDebugPrefix(const std::string& AName)
{
  m_debug_output = AName;
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::execute(const Strategy AStrategy)
{
  auto logger = FrameLogger(m_logEnable);
  logger.log("======================================");
  logger.log(" Starting 2D cross field generation ");
  logger.log("======================================");
  logger.log("Boolean Marks' initialization ");

  initMarks(); 
  logger.log("  -->  DONE");

  //==================================================================
  //STEP 1 - We mark all the nodes, edges, faces classified on
  // boundary entities (points, curves, surfaces)
  //==================================================================
  logger.log("======================================");
  logger.log("Mark boundary cells");
  /* for(auto f_id:m_mesh->faces()){
     std::vector<Node> triangle_verts;
      Face aface = m_mesh->get<Face>(f_id);   
      aface.get<Node>(triangle_verts);
      cout<<"Triangle["<<f_id<<"]= "<<triangle_verts[0].id()<<" "<<triangle_verts[1].id()<<" "<<triangle_verts[2].id()<<endl;  
   }*/
  markBoundaryCells();
  logger.log("  -->  DONE");
    
  //==================================================================
  //STEP 2 - We store all the nodes classified on curves and surfaces
  // in STL vectors
  //=================================================================
  logger.log("======================================");
  logger.log(" Storage of nodes classified on curves ");
  //Mesh::node_iterator it_nodes = m_mesh->nodes_begin();
  m_curve_nodes.clear();
  m_surf_nodes.clear();
 
  for(auto n_id:m_mesh->nodes()){
    Node n = m_mesh->get<Node>(n_id);   

    if (m_mesh->isMarked(n, m_markNodeOnCurv)) 
      m_curve_nodes.push_back(n);
    else if (!m_mesh->isMarked(n, m_markIsolated) &&
	     !m_mesh->isMarked(n, m_markNodeOnPnt) )
      m_surf_nodes.push_back(n);
  }

  logger.log("nodes on curves   : " + std::to_string(m_curve_nodes.size()));
  logger.log("nodes on surfaces : " + std::to_string(m_surf_nodes.size()));
  logger.log("total nb. nodes  : " + std::to_string(m_mesh->getNbNodes()));
  logger.log("    DONE");

  //==================================================================
  //STEP 2 - For nodes on curves, we compute crosses from the 
  // geometric information we have
  //==================================================================
  logger.log("======================================");
  logger.log("Initialization of crosses along curves", FrameLogger::LogEnd::Flush);

  initCrossesOnCurves(); 

  logger.log("    DONE");
  
  //==================================================================
  //STEP 3 - A cross is associated to each node classified on 
  // a geometric point
  //==================================================================
  logger.log("======================================");
  logger.log("Initialization of crosses on points", FrameLogger::LogEnd::Flush);

  initCrossesOnPoints(); 

  logger.log("    DONE");

  if (m_debugEnable)
	writeForDebug();  
  //==================================================================
  // Now, depending on the strategy, the algorithm changes
  //==================================================================
  if(AStrategy==laplace_solve)
    buildFieldViaLaplaceEDP();
  else if(AStrategy== level_sets)
      buildFieldViaLevelSets();
  else if(AStrategy== functional_approach)
          buildFieldViaFunctionalApproach();
  else
    throw GMDSException("Wrong strategy value for 2D cross field generation");
 //==================================================================
  // A final smoothing step can be performed if required
  //==================================================================

  // if(AStrategy!= laplace_solve) {
  //   std::cout << "======================================"<< std::endl;
  //   std::cout << " smoothing" << std::flush;
  //   smoothAll();
  //   std::cout << "    DONE" << std::endl;
  // }
    
  // //==================================================================
  // // DEBUGGING STEP - Tet coloring
  // //==================================================================
  // std::cout << "======================================"<< std::endl;
  // std::cout << " Coloring simplices" << std::endl;
  // colorSimplices();
  // std::cout << "    DONE" << std::endl;

  //==================================================================
  // CLEANING - Boolean marks are cleaned
  //==================================================================
  cleanMarks();
  if (m_debugEnable) {
	writeForDebug();  
	computeReferenceVectorDeviationPerFace();
  }
}

/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::buildFieldViaLaplaceEDP()
{
    std::vector<Face> mesh_faces;
    mesh_faces.resize(m_mesh->getNbFaces());
    int f_index=0; 
    for(auto f_id:m_mesh->faces()){    
        mesh_faces[f_index++]= m_mesh->get<Face>(f_id);
    }
    LaplaceCross2D algo(m_mesh, m_cross_field_2D,
                        m_curve_nodes, m_surf_nodes,
                        mesh_faces);
    algo.execute();
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::buildFieldViaFunctionalApproach()
{/*
    std::vector<Face> mesh_faces;
    mesh_faces.resize(m_mesh->getNbFaces());
    int f_index=0;
    for(auto f_id:m_mesh->faces()){ 
        mesh_faces[f_index++]= m_mesh->get<Face>(f_id);
    }
    FunctionalApproachCross2D algo(this, m_mesh, m_cross_field_2D,
                                   m_curve_nodes, m_surf_nodes,
                                   mesh_faces);
    algo.execute();*/
}
/*---------------------------------------------------------------------------*/

void CrossFieldGeneration2D::buildFieldViaLevelSets()
{/*
  // only nodes on curves are given at the initialization
  // nodes on points have crosses without any meaning
  LevelSetCross2D algo(m_mesh, m_cross_field_2D,
		       m_curve_nodes, m_surf_nodes,
		       m_markFace);
  algo.execute(); */
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::cleanMarks()
{
  m_mesh->unmarkAll<Node>(m_markNodeOnCurv);
  m_mesh->unmarkAll<Node>(m_markNodeOnPnt);
  m_mesh->unmarkAll<Node>(m_markIsolated);
  m_mesh->unmarkAll<Edge>(m_markEdgeOnCurv);
  m_mesh->unmarkAll<Face>(m_markFace);  

  m_mesh->freeMark<Node>(m_markNodeOnCurv);
  m_mesh->freeMark<Node>(m_markNodeOnPnt);
  m_mesh->freeMark<Node>(m_markIsolated);
  m_mesh->freeMark<Edge>(m_markEdgeOnCurv);
  m_mesh->freeMark<Face>(m_markFace); 
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::initMarks()
{     
  m_markNodeOnCurv = m_mesh->newMark<Node>(); 
  m_markNodeOnPnt  = m_mesh->newMark<Node>();
  m_markIsolated   = m_mesh->newMark<Node>();
  m_markEdgeOnCurv = m_mesh->newMark<Edge>();
  m_markFace       = m_mesh->newMark<Face>();  

}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::markBoundaryCells()
{
  BoundaryOperator2D boundaryOp(m_mesh);
  if (!boundaryOp.isValid())
    {
      std::cout << "Invalid model for boundary operations" << std::endl;
      throw GMDSException("Invalid model for boundary operations");
    }
  
  int mark_edge_on_surf = m_mesh->newMark<Edge>();
  int mark_node_on_surf = m_mesh->newMark<Node>();
 
  boundaryOp.markCellOnGeometry(
				m_markEdgeOnCurv, 
				m_markNodeOnCurv, 
				m_markNodeOnPnt,
				m_markIsolated);


  m_mesh->unmarkAll<Node>(mark_node_on_surf);
  m_mesh->unmarkAll<Edge>(mark_edge_on_surf);

  m_mesh->freeMark<Node>(mark_node_on_surf);
  m_mesh->freeMark<Edge>(mark_edge_on_surf);
}
/*---------------------------------------------------------------------------*/
std::vector<Edge> CrossFieldGeneration2D::
getEdgesOnCurve(const Node& ANode) const
{
  std::vector<Edge> edges_on_curve;
  std::vector<Edge> adj_edges = ANode.get<Edge>();
  for (unsigned int i = 0; i < adj_edges.size(); i++)
    {
      Edge ei = adj_edges[i];
      if (m_mesh->isMarked(ei, m_markEdgeOnCurv))
	edges_on_curve.push_back(ei);
    }
  return edges_on_curve;
}
/*---------------------------------------------------------------------------*/
Node CrossFieldGeneration2D::
getNeighboorOn(const Node& ANode, const Edge& AEdge) const
{
  std::vector<Node> nodes = AEdge.get<Node>();
  if (nodes[0].id() == ANode.id())//if (nodes[0].getID() == ANode.getID())
    return nodes[1];

  return nodes[0];
}
/*---------------------------------------------------------------------------*/
void  CrossFieldGeneration2D::initCrossesOnCurves()
{
  //for each node on a geometrical curve, we compute
  //its associated 2D cross
  
  for(auto n_id:m_mesh->nodes()){
      Node current_node = m_mesh->get<Node>(n_id);      
      if (m_mesh->isMarked(current_node, m_markNodeOnCurv) &&
	  !m_mesh->isMarked(current_node, m_markNodeOnPnt)){

	//current_node is on a geometric curve
	//We get adjacent edges that are classified onto
	//a geometric curve. We must get one or two edges
	//at most.
	std::vector<Edge> ridges = getEdgesOnCurve(current_node);
	math::Vector3d newN;
	if (ridges.size() == 1) {
	  Edge current_edge = ridges[0];
	  Node node1 = getNeighboorOn(current_node, current_edge);
	  Node node2 = current_node;

	  //we build the direction vector of the current edge
	  math::Point p1 = node1.point();
	  math::Point p2 = node2.point();
	  math::Vector3d v1 = math::Vector3d(p1, p2);
	  v1.normalize();

	  newN = v1;
	}
	else if (ridges.size() == 2) {
	  //With 2 adajcent edges on the curve, we compute average values

	  Edge edge1 = ridges[0];
	  Edge edge2 = ridges[1];

	  Node node1 = getNeighboorOn(current_node, edge1);
	  Node node2 = getNeighboorOn(current_node, edge2);

	  //we build the average between adjacent edges
	  math::Point  p1 = node1.point();
	  math::Point  p  = current_node.point();
	  math::Point  p2 = node2.point();
	  math::Vector3d v1 = math::Vector3d(p1, p);
	  math::Vector3d v2 = math::Vector3d(p, p2);
	  v1.normalize();
	  v2.normalize();

	  newN = v1+v2;
	  newN.normalize();
	}
	else{
	  std::cout << "Nb ridges for an edge adjacent to node "
		    << current_node.id() << ": " << ridges.size() << std::endl;
	  throw GMDSException("A ridge node has an illegal number of edges.");
	}
	std::vector<Face> adj_faces;
	current_node.get<Face>(adj_faces);
	Face adj_f = adj_faces[0];
	math::Vector3d normal  = adj_f.normal();
	
	math::Vector3d v1 = newN;

	//Computation of v3 from v1 and v2.
	math::Vector3d v2 = v1.cross(normal);
	math::Cross2D c(v1,v2);

	(*m_cross_field_2D)[current_node.id()] = c;
      }//if (m_mesh->isMarked(current_node, m_markNodeOnCurv))

    }//for (; !it_node.isDone(); it_node.next())
}
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::initCrossesOnPoints()
{
   for(auto n_id:m_mesh->nodes()){
    Node current_node = m_mesh->get<Node>(n_id);
      if (m_mesh->isMarked(current_node, m_markNodeOnPnt))
	{
	  //we initialize the value at points, but we do not use it after.
	  (*m_cross_field_2D)[current_node.id()] = math::Cross2D();
	}

    }
}
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::smoothAll()
{
  int mark_smooth = m_mesh->newMark<Node>();
  for(auto n_id:m_mesh->nodes()){
    Node n = m_mesh->get<Node>(n_id);
      if (!m_mesh->isMarked(n, m_markNodeOnCurv) &&
	  !m_mesh->isMarked(n, m_markNodeOnPnt))
	m_mesh->mark(n, mark_smooth);

    }
  smooth(mark_smooth);

  m_mesh->unmarkAll<Node>(mark_smooth);
  m_mesh->freeMark<Node>(mark_smooth);



} 
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::smooth(const int AMark)
{

  //	SmoothingHLBFGS smoother(m_mesh, m_cross_field_2D, m_surf_normal);
  /*	FrameFieldLaplacianSmoothing smoother(m_mesh, m_cross_field_2D, m_surf_normal);
	smoother.initBoundaryMarks(m_markNodeOnPnt, m_markNodeOnCurv, m_markNodeOnSurf);
	smoother.selectNodes(AMark);
	smoother.execute();
  */
}
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::colorSimplices()
{
  // Variable<int>* var_sing = 0;
  // try{
  //   var_sing = m_mesh->getVariable<int>(GMDS_REGION, "sing_tet");
  // }
  // catch (GMDSException& e){
  //   var_sing = m_mesh->newVariable<int>(GMDS_REGION, "sing_tet");
  // }
		    
  // IGMesh::region_iterator it = m_mesh->regions_begin();
		      
  // int nbOfCluster = 0;
  // int nbColoredTet = 0;
  // //=========================================================================
  // // INTERN SKELETON CREATION
  // //=========================================================================
  // // FIRST LOOP ON REGIONS TO GET ALL THE 3-SING. TETS
  // //=========================================================================
  // for (; !it.isDone(); it.next()){
  //   Region current_region = it.value();
  //   std::vector<Node> nodes = current_region.get<Node>();
  //   bool onPnt = false;
  //   for (unsigned int i_node = 0; i_node < nodes.size(); i_node++)
  //     {
  // 	Node ni = nodes[i_node];
  // 	if (m_mesh->isMarked(ni, m_markNodeOnPnt))
  // 	  onPnt = true;

  //     }
  //   if (onPnt)
  //     {
  // 	(*var_sing)[current_region.getID()] = 0;
  //     }
  //   else
  //     {
  // 	std::vector<TCellID> nodeIDs = current_region.getIDs<Node>();
  // 	int ID1 = nodeIDs[0];
  // 	int ID2 = nodeIDs[1];
  // 	int ID3 = nodeIDs[2]; 
  // 	int ID4 = nodeIDs[3];
						
  // 	int singTypeTmp = math::Quaternion::testSingularity((*m_cross_field)[ID1],
  // 							    (*m_cross_field)[ID2],
  // 							    (*m_cross_field)[ID3],
  // 							    (*m_cross_field)[ID4]);
  // 	if (singTypeTmp != 0)
  // 	  nbColoredTet++;
						    
  // 	(*var_sing)[current_region.getID()] = singTypeTmp;
  //     }
  // }
  // std::cout << "Nb. colored tetrahedra: " << nbColoredTet << std::endl;
}

/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::writeForDebug(const std::string AFileName)
{
  static int nb_file = 0;
  Variable<double>* var_angle = 0;
  Variable<math::Vector3d>* var_X = 0;
  Variable<math::Vector3d>* var_Y = 0;
  Variable<math::Vector3d>* var_XM = 0;
  Variable<math::Vector3d>* var_YM = 0;
  try{
    var_angle = m_mesh->getVariable<double,GMDS_NODE>("angle");
    var_X     = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("cross_X");
    var_Y     = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("VY");
    var_XM    = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("VXM");
    var_YM    = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("VYM");
  }
  catch (GMDSException& e){
    var_angle = m_mesh->newVariable<double,GMDS_NODE>("angle");
    var_X     = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("cross_X");
    var_Y     = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("VY");
    var_XM    = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("VXM");
    var_YM    = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("VYM");
  }
 
 for(auto n_id:m_mesh->nodes()){
    math::Cross2D ci = (*m_cross_field_2D)[n_id];
    (*var_angle)[n_id] = ci.referenceAngle();
    std::vector<math::Vector3d> current_vectors = ci.componentVectors(); 
    (*var_X)[n_id] = current_vectors[0];
    (*var_Y)[n_id] = current_vectors[1];
    (*var_XM)[n_id] = current_vectors[2];
    (*var_YM)[n_id] = current_vectors[3];
 }
  gmds::IGMeshIOService meshIoServ(m_mesh);
  gmds::VTKWriter writer(&meshIoServ);
    writer.setCellOptions(gmds::N|gmds::F);
    writer.setDataOptions(gmds::N|gmds::F);
  //VTKWriter<Mesh> writer(*m_mesh);
  if (AFileName != "")
    { 
      std::stringstream file_name;
      file_name <<m_debug_output<<"_"<<AFileName;     
      writer.write(file_name.str()); // writer.write(file_name.str(), DIM3 | F | N);
    }
  else
    {
      std::stringstream file_name;
      file_name <<m_debug_output<<"FFG_2D_Debug_" << nb_file<<".vtk";;
      std::cout<<file_name.str()<<std::endl;
      writer.write(file_name.str());// writer.write(file_name.str(), DIM3 | F | N);
      nb_file++;
    }

  double x_min = 100000;
  double y_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
  
  for(auto n_id:m_mesh->nodes()){
    Node n = m_mesh->get<Node>(n_id);    
      math::Point p = n.point();
      if (p.X() < x_min)
	x_min = p.X();
      if (p.X() > x_max)
	x_max = p.X();

      if (p.Y() < y_min)
	y_min = p.Y();
      if (p.Y() > y_max)
	y_max = p.Y();


    }
  double dist_x = x_max - x_min;
  double dist_y = y_max - y_min;

  double cube_size = 0;
  if (dist_x <= dist_y ){
    cube_size = dist_x;
  }
  else {
    cube_size = dist_y;
  }

  cube_size /= 20;

  MeshModel model_cube(DIM3 | F | N | F2N);
  Mesh mesh_cube(model_cube);
  
  for(auto n_id:m_mesh->nodes()){
    Node n = m_mesh->get<Node>(n_id);
      math::Point center = n.point();
      //      if (m_mesh->isMarked(n, m_mark_alive))
{
	math::Cross2D current_cross = (*m_cross_field_2D)[n.id()];

	std::vector<math::Vector3d> current_vectors = current_cross.componentVectors();
	math::Vector3d vx = current_vectors[0];
	math::Vector3d vy = current_vectors[1];
	math::Point p1 = center + (vx + vy )*cube_size;
	Node n1 = mesh_cube.newNode(p1);
	math::Point p2 = center + (vx - vy)*cube_size;
	Node n2 = mesh_cube.newNode(p2);
	math::Point p3 = center + (vx + vy).opp()*cube_size;
	Node n3 = mesh_cube.newNode(p3);
	math::Point p4 = center + (vy - vx)*cube_size;
	Node n4 = mesh_cube.newNode(p4);
	mesh_cube.newQuad(n1, n2, n3, n4);
      }
    }
  
  gmds::IGMeshIOService meshIoServCube(&mesh_cube);
  gmds::VTKWriter writer_cube(&meshIoServCube);
  writer_cube.setCellOptions(N|F);
  writer_cube.setDataOptions(N|F);
  
  std::stringstream file_name_cube;
  file_name_cube<<m_debug_output <<"FFG_2D_Debug_Cube_" << nb_file<<".vtk";
 
  writer_cube.write(file_name_cube.str());


}


/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::computeReferenceVectorDeviationPerFace(){
    
     Variable<double>* var_alpha = 0;
     var_alpha     = m_mesh->newVariable<double,GMDS_FACE>("angleDeviation");
	Variable<double>* var_alpha_max = 0;
     var_alpha_max     = m_mesh->newVariable<double,GMDS_FACE>("maxAngleDeviation");
     
    for(auto f_id:m_mesh->faces()){
        vector<gmds::Node> current_nodes = m_mesh->get<Face>(f_id).get<Node>();
             
        math::Cross2D cross_0 = (*m_cross_field_2D)[current_nodes[0].id()];
        math::Cross2D cross_1 = (*m_cross_field_2D)[current_nodes[1].id()];
        math::Cross2D cross_2 = (*m_cross_field_2D)[current_nodes[2].id()];
        
        math::Vector3d RV0 = cross_0.referenceVector();
        math::Vector3d RV1 = cross_1.referenceVector();
        math::Vector3d RV2 = cross_2.referenceVector();
        
	   (*var_alpha)[f_id] = RV0.angle(RV1) + RV1.angle(RV2) + RV2.angle(RV0);
        (*var_alpha)[f_id] = (*var_alpha)[f_id]/3.0;
	   (*var_alpha_max)[f_id] = std::max(std::max(RV0.angle(RV1), RV1.angle(RV2)), RV2.angle(RV0));
    }
    
    
  gmds::IGMeshIOService meshIoServCube(m_mesh);
  gmds::VTKWriter writer_cube(&meshIoServCube);
  writer_cube.setCellOptions(N|F);
  writer_cube.setDataOptions(N|F);
  
  std::stringstream file_name;
  file_name<<m_debug_output<<"-angleDeviation.vtk";
 
  writer_cube.write(file_name.str());
}
