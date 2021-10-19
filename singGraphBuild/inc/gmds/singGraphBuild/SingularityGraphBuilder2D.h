/*----------------------------------------------------------------------------*/
/*
 * SingularityGraphBuilder2D.h
 *
 *  Created on: April 10 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef SINGULARITYGRAPHBUILDER_2D_H_
#define SINGULARITYGRAPHBUILDER_2D_H_
/*----------------------------------------------------------------------------*/
#include <cstdlib>
#include <iostream>
#include <string>
#include <list>
#include <map>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Cross2D.h>
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingularityGraph.h>
#include <gmds/singGraphBuild/Tools.h>
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
/*----------------------------------------------------------------------------*/
//#include <Tools.h>
/*----------------------------------------------------------------------------*/
/** \brief Class providing an algorithm to build a 2D singularity graph from
 *         a triangular mesh and a 2D cross field defined on this mesh
 */
class LIB_GMDS_SINGGRAPHBUILD_API SingularityGraphBuilder2D
{
public:
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Strategies to compute the singularity graph
	* \param original - the original version will be implemented: start from a slot and advance from triangle to triangle using the Heun method, constructing thus the singularity line;
	* I. if the singularity line is within the confusing ball of a different singularity slot, we connect the two slots (if no previous "better"(smaller deviation) connection for that slot has already been detected) and add the line
	*       to the singularity graph
	* II. if the singularity line reached the boundary, we add it to the final singularity graph 
	* \param simultaneousStartHeun - the simultaneous strategy is implemented; we simultaneously depart from all slots of all the singularities (advancing from triangle to triangle using Heun's method)
	* we therfore construct singularity lines from all slots; after advancing each such line we test whether we are within a certain distance (thresholdStreamLineDist) from the ending point of a different singularity line, we connect the two.  
	* \param simultaneousStartRK4 - same sa the previous strategy, except we advance using RK4 - this has the advantage of being more accurate; we no longer advance from triangle to neighbouring triangle (like for previous strategies), 
	* but rather we advance by a certain geometric distance (step), having therefore the same "speed" for all simultaneuos departed lines.
	* \param shortestPaths - this strategy is a discrete one: we use a modiffied Dijkstra's algorithm to find the "shortest paths" (as having the least deviation from the cross field);
	* for each singularity point we advance using RK4 until the end points of all the slots of that singularity are situated in triangles that are not neighbouring (nor by edge, nor by node)
	*  afterwards, for each slot of each singularity we try to find the shortest paths from it to all other slots of all other singularities and we assign them a weight (inverse proportional to the allignement of the shortest path line to the field).
	* NOTE: for now, we do not allow that a slot of one singularity can be connected to a different slot of the same singularity; however this could be envisaged in certain configurations
	* Having detected all the possible shortest path lines from all slots to all other slots we assemble them into a system with constraints that will indicate the best possible configuration;
	* the system that we must solve has binary variables: the existence (or not) in the final graph of each shortest path line and the following constraints: if there are a number of N slots, the total number of variables is N*(N+1) (accounting also for shortest path lines to boundary); 
	* NOTE: indeed we could have less variables in the system (slot x towards slot x shouldn't exist and also slot x towards slot y should be the same as slot y towards slot x); however, in our test, it seems the solve of the system is quite fast, the assembling is the one that takes the longest time (the actual Dijkstra's step)
	* I. The objective function f =  sum (var * weight[variable]).
	* II. The constraints: 
	* a) each slot has exactly one shortest path line departing from it.
	* b) two (or more)  shortest path lines that pass through the same triangle and follow the same cross component (for simplicity in 2D those are vertical or horizontal) are not allowed to co-exist in the final graph.
	* - these are detected and stored into "IllegalCross".
	* c) each slot can be connected only to slots of different singularities (see the note above).
	* d) if at the initial step (advancing with RK4 from a singilarity point in the direction of all its slots until reaching non-neighbouring triangles) we reach the boundary or a different slot, the variable for this connection 
	*will be set to true and will be fixed.      
	* \param testPrescribedSing - performs the shortestPaths strategy, but starting from singularity triangles that are manually introduced by the user (their ID) - this is present for testing purposes
	*/
	typedef enum {
	             original,        
	             simultaneousStartHeun,
	             simultaneousStartRK4,
	             shortestPaths,
	             testPrescribedSing
	} Strategy;
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Constructor.
	* \param AMesh the mesh where we work on
	* \param AField the cross field associated to AMesh
	* \param ABuildGeomSing flag to build the geometric singularity points.
	*/
	SingularityGraphBuilder2D(gmds::Mesh* AMesh,
	                          gmds::Variable<gmds::math::Cross2D>* AField,
	                          const bool ABuildGeomSing = true);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Execution of the algorithm with the chosen AStrategy
	*/
	void execute(const Strategy AStrategy, unsigned int number_of_control_points);
    /*----------------------------------------------------------------------------------------------------*/
    /** \brief gives acces to the build graph
    */
	const SingularityGraph& graph() const {return m_graph;}
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Function to give the directory where we want to put the output
	*		   files (vtk files)
	*/
	void setDebugPrefix(const std::string& ADirName) {
		m_output_directory_name = ADirName;
	}
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Initialization of the boolean marks used for the algorithm
	* \param AMarkNodePnt mark for nodes classified on points
	* \param AMarkNodeCrv mark for nodes classified on curves
	* \param AMarkEdgeCrv mark for edges classified on curves
	*/
	void initMarks(const int AMarkNodePnt, const int AMarkNodeCrv,
	               const int AMarkEdgeCrv);

	/*----------------------------------------------------------------------------------------------------*/
    
protected:
	void colorFaces(const int markF, const int markE);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief  Detect the triangles that contain singularity points
	*/
	void detectSingularTriangles();
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief  For a a face AFace, which is singular, this method computes the
	* geometric location of the singularity point in AFace
	* \param[in] AFace    a singular face
	* \param[out] APosSing the singularity point location
	*/
	void computeSingPointInfo(const gmds::Face&        AFace,
	                          gmds::math::Point& APosSing);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief  For a a face AFace, which is singular, this method creates the
	* the singularity point and the slots (cont - used for visualization purposes)
	*/
	void createSingPointAndSlots(const gmds::Face& AFace);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Add the boundary nodes and edges into the singularity graph.
	* \param[in] artificialSingPointsCreated  a vector containing the artifial singularity points 
	* that have been created (for each "hole" in the mesh we create one such point)
	*/
	void addGeometryToSingularityGraph(vector<CurveSingularityPoint*> &artificialSingPointsCreated);
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Build geometric slots. As a consequence, it means the resulting
	* domain will be split and  made of 4-sided patches only.
	*/
	void buildGeometricSlots();
	
	/*----------------------------------------------------------------------------------------------------*/
	void createLineFrom(VertexSingularityPoint* AFrom, 
	                    const double            AAngle,
	                    const int               ANbLines);
	/*----------------------------------------------------------------------------------------------------*/
	
	void writeOutput(const std::string& AFileName);
	/*----------------------------------------------------------------------------------------------------*/
	
	void writeOutputSingle(const std::string& AFileName);
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having a block structure dictated by the resulted 
	*graph; This mesh will have as vertices the singularity points (and possibly geometric points), 
	* as well as the intersection of the resulted singularity lines with one another and with the boundary. 
	* Normally this function creates a collection of patch elements that are rectangular (only 4 verices and 4 edges). 
	* However, if the user wants to visualize the patches having as boundary the original detected singularity lines, 
	* it can set "bool curvePatches = true;" inside the function.
	*/
	void writeOutputPatches(const std::string& AFileName);
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having as vertices some of the original mesh's vertices \param ANodes and zero-are triangles
	* therefore writes a vtk file containing lines
	*/
	void writeTestMeshVerts(vector<gmds::Node>& ANodes, std::string& AFileName);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having as vertices \param APoints
	* and 'zero-area' triangles
	*/
	void writeTestPoints(vector<gmds::math::Point>& APoints, std::string& AFileName);
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having as triangles \param ATriangles
	*/
	void writeTestMeshTriangles(vector<gmds::Face>& ATriangles, std::string& AFileName);    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having as triangles a subset of the original mesh's triangles (\param ATrianglesIds), given as ids
	*/
	void writeTestMeshTrianglesIds(vector<gmds::TCellID>& ATrianglesIds, std::string& AFileName); 
	/*----------------------------------------------------------------------------------------------------*/
	/* write a vtk file containing the confusing balls*/
	void writeConfusingBalls();
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creation of singularity points
	*/
	//  void createSingularityPoints();
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creation of singularity lines - original strategy
	*/
	virtual void createSingularityLines() = 0;
	/*----------------------------------------------------------------------------------------------------*/
    
	/** \brief Remove a singularity line
	* \param ALine  the line that must be removed from the graph
	*/
	void removeSingularityLine(SingularityLine* ALine);

	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Compute the singularity line starting from point \param AFromSingPnt in the
	*         direction of \param AFromSlot as long as the accumulated distance 
	*         is smaller than the the \param searchStep or until it hasn't reached the 
	* confusing ball of a different singularity or until the line reaches the boundary (using Heun)
	*
	* \param[in] AFromPnt      the singularity point we start from
	* \param[in] AFromSlot     the associated slot we really consider
	* \param[out] AToSingPnt    the singularity point we arrive at, if it exists
	* \param[out] AToSlot       the associated slot
	* \param[out] AToPnt        the location we arrive at the end
	* \param[out] AToDir        the last direction to arrive at the end
	* \param[out] APoints       a discretization of the stream line
	* \param[out] ATriangles    the list of traversed triangles
	* \param[out] AToCellDim    the dimension of the cell we finish on (0,1,2)
	* \param[out] AToCellID     the id of the cell we finish on
	* \param[out] streamlineDeviation     the deviation of the computed streamline     
	* \param[out] accumulatedDistancePerSlotLine     the geometric distance of the line
	* \param[out] searchStep    step to grow the line (geometric length)
	* \param[out] AEndOnBnd     indicates if we finish on the boundary (true)
	* \param[out] AToSlotIsFree indicates if we finish onto a free slot (true)
	*/
	/*----------------------------------------------------------------------------------------------------*/
	void growLine(SingularityPoint*               AFromSingPnt,
	              SingularityPoint::Slot*         AFromSlot,
	              SingularityPoint*&              AToSingPnt,
	              SingularityPoint::Slot*&        AToSlot,
	              gmds::math::Point&              AFromPnt,
	              gmds::math::Vector3d&           AToDir,
	              std::vector<gmds::math::Point>& APoints,
	              std::vector<gmds::TCellID>&     ATriangles,
	              int&                            AToCellDim,
	              gmds::TCellID&                  AToCellID,
	              double&                         streamlineDeviation,                  
	              double&                         accumulatedDistancePerSlotLine,
	              double&                         searchStep,
	              bool&                           AEndOnBnd,
	              bool&                           AToSlotIsFree);
    
	/*----------------------------------------------------------------------------------------------------*/
   
	/** \brief Compute the singularity line starting from point \param AFromSingPnt in the
	*         direction of \param AFromSlot by \param stepSize as long as the accumulated distance 
	*         is smaller than the the \param searchStep or until it hasn't reached the 
	* confusing ball of a different singularity or until the line reaches the boundary (using RK4)
	*
	* \param[in] AFromPnt      the singularity point we start from
	* \param[in] AFromSlot     the associated slot we really consider
	* \param[out] AToSingPnt    the singularity point we arrive at, if it exists
	* \param[out] AToSlot       the associated slot
	* \param[out] AToPnt        the location we arrive at the end
	* \param[out] AToDir        the last direction to arrive at the end
	* \param[out] APoints       a discretization of the stream line
	* \param[out] ATriangles    the list of traversed triangles
	* \param[out] AToCellDim    the dimension of the cell we finish on (0 -Node, 1- Edge, 2 - Triangle)
	* \param[out] AToCellID     the id of the cell we finish on
	* \param[out] streamlineDeviation     the deviation of the computed streamline     
	* \param[out] accumulatedDistancePerSlotLine     the geometric distance of the line
	* \param[out] searchStep    step to grow the line 
	* \param[out] AEndOnBnd     indicates if we finish on the boundary (true)
	* \param[out] AToSlotIsFree indicates if we finish onto a free slot (true)
	* \param[out] find_end     indicates if we have reached the end (either boundary, a different slot or if we cannot continue) (true)
	*/
	void growLineRK4(const SingularityPoint::Slot*   AFromSlot,
	                 SingularityPoint::Slot*&        AToSlot,
	                 const gmds::math::Point&        AFromPnt,
	                 gmds::math::Vector3d&           AToDir,
	                 std::vector<gmds::math::Point>& APoints,
	                 std::vector<gmds::TCellID>&     ATriangles,
	                 int&                            AToCellDim,
	                 gmds::TCellID&                  AToCellID,
	                 double&                         streamlineDeviation,   
	                 const double&                   stepSize,
	                 bool&                           AEndOnBnd,
	                 bool&                           AToSlotIsFree,
	                 bool&                           find_end);
    
	/*----------------------------------------------------------------------------------------------------*/    
    
	/** \brief Initialize the confusing ball for singularity \param APnt
	*
	* \param APnt the singularity point we work on
	*/
	void initConfusingBalls(SingularityPoint* APnt);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Redefine (usually grow) the radius of the confusing ball for singularity \param APnt
	* by \param increaseRadiusScale and store the previous radius \param previousRad  
	* \param APnt the singularity point we work on
	*  \param modifiedFaces the list of modiffied faces (those triangles that are now situated within the new radius)
	*  \param previousRad the previous radius
	*  \param increaseRadiusScale the amount by which we increase the previous radius
	*/
    
	void redefineOneConfusingBall(SingularityPoint* APnt, vector<gmds::TCellID>& modifiedFaces, double &previousRad, const double increaseRadiusScale = 1.5);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creates a geometric singularity point.
	*
	* \param [in]  AInPnt   the point where we are located
	* \param [in]  AInVec   the geometric direction we arrive
	* \param [in]  ACellDim the dimension of the cell we are located on
	* \param [in]  ACellID  the id of the cell we are located on
	* \param [out] APnt the created singularity point
	* \param [out] ASlot the associated slot directed towards AInVec
	*/
	void createGeometricSingularityPoint(const gmds::math::Point&    AInPnt,
	                                     const gmds::math::Vector3d& AInVec,
	                                     const int                   ACellDim,
	                                     const gmds::TCellID         ACellID,
	                                     SingularityPoint*&          APnt,
	                                     SingularityPoint::Slot*&    ASlot);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Returs the singularity lines that go througt the face \param AFace
	*
	* \param AFace a mesh face
	* \return the singularity lines going through AFace
	*/
	std::vector<SurfaceSingularityLine*> getSingularityLinesIn(const gmds::Face& AFace);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief This operation detects line intersection in steady areas and
	*         creates singularity points
	*/
	virtual void detectLineIntersections();
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief This operation creates intersections between 2 lines in the
	*         vicinity of a mesh face
	*/
	void createLineIntersection(SurfaceSingularityLine*         ALine1,
	                            SurfaceSingularityLine*         ALine2,
	                            gmds::Face&                     AFace,
	                            std::vector<gmds::math::Point>& AAddedPoints);
	/*----------------------------------------------------------------------------------------------------*/
	void createLineIntersection(std::vector<SurfaceSingularityLine*>& ALines,
	                            gmds::Face&                           AFace,
	                            std::vector<gmds::math::Point>&       AAddedPoints);    
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Computes the barycenters for all triangles in the mesh, 
	* as well as the crosses associated to them
	*/
	void constructCenterTriangleCrosses();
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Computes the barycenter (\param centerTri) for triangle \param ATriangle, 
	*as well as the cross associated to it (\param centerTriCross) */
	void constructOneCenterTriangleCross(gmds::Face& ATriangle, 
	                                     gmds::math::Point &centerTri, 
	                                     gmds::math::Cross2D& centerTriCross);
    
	/*----------------------------------------------------------------------------------------------------*/
     
	/** \brief Constructs a new mesh (\param newLocalMesh) that has as Nodes all Nodes of 
	* \param trianglesToRemesh plus their centers and the resulted new triangles 
	* (each original triangle splits into 3 triangles)
	* \param newLocalMesh the newly constructed mesh    
	* \param newLocalMesh_id_to_mesh_id_node in this vector at the position corresponding the Node id of the \param newLocalMesh (except new introduced vertices) it is stored the id of the corresponding vertex in the original mesh 
	* \param local_cross_field_2D cross field of the \param newLocalMesh
	* \param trianglesToRemesh the ids of the triangles of the original mesh that we want to remesh
	*/
	void remeshTriangles(gmds::Mesh*                          newLocalMesh,  
	                     vector<gmds::TCellID>&               newLocalMesh_id_to_mesh_id_node, 
	                     gmds::Variable<gmds::math::Cross2D>* local_cross_field_2D, 
	                     vector<gmds::TCellID>&               trianglesToRemesh);
     
     
	/*----------------------------------------------------------------------------------------------------*/
     
	/** \brief Constructs a new mesh (\param newLocalMesh) as a copy of the original mesh, 
	 * except some triangles (\param trianglesToRemesh)
	*
	* \param newLocalMesh the newly constructed mesh
	* \param newMesh_cross_field_2D cross field of the \param newLocalMesh
	* \param trianglesToRemesh the ids of the triangles of the original mesh that we want to remesh
	* \param remeshedTriangles boolean value indicating whether a triangle has been or should be remeshed
	* \param newTriangleCenters the centers of all triangles of the \param newLocalMesh
	* \param newBdryEdgeNormals the normals for all new boundary edges
	* \param newBdryNodeNormals the normals for all new boundary nodes
	* \param isCurveEdge boolean value indicating whether an edge is a boundary edge
	* \param isCurveNode boolean value indicating whether a node is a boundary node
	* \param newLocalMesh_id_to_mesh_id_node in this vector at the position corresponding the Node id of the \param newLocalMesh (except newly introduced vertices) it is stored the id of the corresponding vertex in the original mesh 
	* \param mesh_id_to_newLocalMesh_id_node in this vector at the position corresponding the Node id of original mesh it is stored the id of the corresponding vertex in the \param newLocalMesh
	* \param newLocalMesh_id_to_mesh_id_face in this vector at the position corresponding the Face id of the \param newLocalMesh (except newly introduced Face) it is stored the id of the corresponding Face in the original mesh 
	* \param mesh_id_to_newLocalMesh_id_face in this vector at the position corresponding the Face id of original mesh it is stored the id of the corresponding Face in the \param newLocalMesh
	*/
	void remeshTrianglesNewMesh(gmds::Mesh*                          newLocalMesh,   
	                            gmds::Variable<gmds::math::Cross2D>* newMesh_cross_field_2D, 
	                            vector<gmds::TCellID>&               trianglesToRemesh,
	                            vector<bool>&                        remeshedTriangles,
	                            vector<gmds::math::Point>&           newTriangleCenters, 
	                            vector<gmds::math::Cross2D>&         newTriangleCenterCrosses, 
	                            vector<gmds::math::Vector3d>&        newBdryEdgeNormals,
	                            vector<gmds::math::Vector3d>&        newBdryNodeNormals,
	                            vector<bool>&                        isCurveEdge,
	                            vector<bool>&                        isCurveNode,
	                            vector<gmds::TCellID>&               newLocalMesh_id_to_mesh_id_node,
	                            vector<gmds::TCellID>&               mesh_id_to_newLocalMesh_id_node,
	                            vector<gmds::TCellID>&               newLocalMesh_id_to_mesh_id_face,
	                            vector<gmds::TCellID>&               mesh_id_to_newLocalMesh_id_face);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Connect two singularity lines if they are within the given threshold distance; 
	* if(connectByField) connect them if by folowing the frame field we will not change direction (closest component vector along the connection path will have the same direction 0-2 or 1-3 from start to end)
	*
	* \param[in] dir_slot_i     the last direction of singline1 before deciding to connect
	* \param[in] dir_slot_j     the last direction of singline2 before deciding to connect
	* \param[in] connection_dir the direction of the straight line connecting last points of singline1 and singline2
	* \param[in] start_pnt1     the last point of singline1 before deciding to connect
	* \param[in] start_pnt2     the last point of singline2 before deciding to connect   
	* \param[in] start_cell_id1     the last cell id of singline1 before deciding to connect 
	* \param[in] start_cell_id2     the last cell id of singline2 before deciding to connect
	* \param[in] start_cell_dim1     the last cell dimension of singline1 before deciding to connect 
	* \param[in] start_cell_dim2     the last cell dimension of singline2 before deciding to connect
	* \param[out] APoints       a discretization of the stream line
	* \param[out] ATriangles    the list of traversed triangles    
	* \param[out] streamlineDeviation     the deviation of the computed streamline
	*/
     
	void connectSingularityLines(gmds::math::Vector3d            dir_slot_i,  
	                             gmds::math::Vector3d            dir_slot_j,
	                             gmds::math::Vector3d            connection_dir,
	                             gmds::math::Point               start_pnt1,// start_pnt = line_discretization[0]; //starting point
	                             gmds::math::Point               start_pnt2, 
	                             gmds::TCellID                   start_cell_id1,
	                             gmds::TCellID                   start_cell_id2,
	                             int                             start_cell_dim1, 
	                             int                             start_cell_dim2,
	                             std::vector<gmds::math::Point>& APoints,
	                             std::vector<gmds::TCellID>&     ATriangles,                 
	                             double&                         streamlineDeviation);
	
	/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/** \brief get the shortest path between a singular face and all other singular faces by 'walking from face center to face center' in the new local mesh
	*
	* \param[in] source    the singular face from which we start
	* \param[in] targets   all the other singular face
	* \param[out] min_distance    min_distance[i] - minimum distance between source and face "i" - 
	* 						first - total distance
	* 						second - number of crossed faces to get to "i" from source
	* 						this distance is a combination between the matching numbers between 2 faces and the 
	* deviation of the cross field (across the centers of all visited faces)
	* \param[out] previous     the id of the previous visited face
	* \param[out] found        the first singularity face we arrive into
	* \param[in] face2FaceNeighboursByVerts     the neighbours of a face (considered as the adjacent faces to the 3 vertices of the current face)
	* \param[in] face2FaceTransport     the matching number between a face and its neighbours (same order as face2FaceNeighbours) 
	* \param[in] face2FaceDeviation     the deviation between a face and its neighbours (same order as face2FaceNeighbours) (calculated as the deviation of the cross field in between triangle centers)
	* \param[in] isBdryFace       boolean value indicating if the face is on the boundary  
	* \param[in] prevDir          previous direction
	* \param[in] prevDirCross     previous direction and previous cross (dictated by the previous triangle in which we were)
	* \param[in] targetDirCross   target (opposite) direction and cross (dictated by the triangle in which we want to arrive)
	* \param[in] targetBdry boolean value indicating if our target is the boundary
	* \param[in] startPoint last added point (using RK4) for the begining slot; it is located inside face \param[in] source
	* \param[in] endPoint last added point (using RK4) for the enfing slot; it is located inside face \param[in] targets[0];
	* \param[out] endPoint last added point on the boundary
	* \param[in] maxDist  minimum weight for penalizing a path which switched the cross component direction 
	*/
	void getShortestPathBtwFacesOptimizedNewLocalMesh(gmds::TCellID&                                                              source,
	                                                  vector<unsigned int>&                                                       targets,
	                                                  vector<pair<double, unsigned int>>&                                         min_distance,
	                                                  vector<int>&                                                                previous, 
	                                                  int&                                                                        found,
	                                                  vector<vector<double>>&                                                     face2FaceMatch, 
	                                                  vector<vector<unsigned int>>&                                               face2FaceDeviation,
	                                                  pair<gmds::math::Vector3d, gmds::math::Vector3d>&                           prevDirCross,
	                                                  vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>>&                   targetDirCross,
	                                                  bool&                                                                       targetBdry,
	                                                  vector<bool>&                                                               is_bdry_face,
	                                                  gmds::math::Point&                                                          startPoint,
	                                                  gmds::math::Point&                                                          endPoint,
	                                                  double&                                                                     maxDist,
	                                                  int&                                                                        to_cell_dim, 
	                                                  gmds::TCellID&                                                              to_cell_id,
	                                                  gmds::TCellID&                                                              lastVisTri,
	                                                  unsigned int&                                                               contSource,
	                                                  unsigned int&                                                               colNo,
	                                                  unsigned int&                                                               totalNumberOfSlots,
	                                                  vector<unsigned int>&                                                       faceNo2Cont,
	                                                  vector<vector<vector<gmds::math::Point>>>&                                  line_discretization,
	                                                  vector<vector<vector<gmds::TCellID>>>&                                      line_triangles,
	                                                  vector<vector<pair<vector<gmds::math::Point>, vector<gmds::math::Point>>>>& pointPaths,
	                                                  vector<vector<vector<gmds::TCellID>>>&                                      finalPaths,
	                                                  vector<vector<vector<gmds::math::Point>>>&                                  finalCenterPoints,
	                                                  vector<vector<vector<gmds::TCellID>>>&                                      traversedSPTris,
	                                                  vector<vector<double>>&                                                     distances,
	                                                  vector<pair<unsigned int, unsigned int>>&                                   IllegalCross,
	                                                  vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedFace,
	                                                  vector<pair<unsigned int, unsigned int>>&                                   contSourceToSingularity,
	                                                  vector<bool>&                                                               singOrGeomFaces,
	                                                  gmds::Mesh*                                                                 newLocalMesh,
	                                                  vector<gmds::TCellID>&                                                      newLocalMesh_id_to_mesh_id_node,
	                                                  vector<gmds::TCellID>&                                                      mesh_id_to_newLocalMesh_id_node,
	                                                  vector<gmds::TCellID>&                                                      newLocalMesh_id_to_mesh_id_face,
	                                                  vector<gmds::TCellID>&                                                      mesh_id_to_newLocalMesh_id_face,
	                                                  gmds::Variable<gmds::math::Cross2D>*                                        newMesh_cross_field_2D,
	                                                  vector<bool>&                                                               trianglesToRemeshBool,
	                                                  vector<gmds::math::Vector3d>&                                               newBdryEdgeNormals,
	                                                  vector<gmds::math::Vector3d>&                                               newBdryNodeNormals,
	                                                  vector<bool>&                                                               isCurveEdge,
	                                                  vector<bool>&                                                               isCurveNode,
	                                                  vector<gmds::math::Point>&                                                  newTriangleCenters, 
	                                                  vector<gmds::math::Cross2D>&                                                newTriangleCenterCrosses);
	/*----------------------------------------------------------------------------------------------------*/
    
	/** \brief get the shortest path between a singular face and all other singular faces by walking through a 
	* mesh whose nodes are the nodes of the original triangles + the centers of the original triangles;
	* the ids for the new mesh are as follows: the first ids correspond to the original nodes and the following 
	* are in the interval[m_mesh->getNbNodes(), m_mesh->getNbNodes() + m_mesh->getNbFaces()], in the order of the face ids
	*
	* \param[in] numberOfNewNodes     total number of nodes in the new mesh (m_mesh->getNbFaces() + m_mesh->getNbNodes())
	* \param[in] source    the singular nodes from which we start (the center of a slot face)
	* \param[in] targets   all the other singular nodes (the centers of all the other slot faces - except the slots of the source slot singularity)
	* \param[out] min_distance    min_distance[i] - minimum distance between source and node "i" - 
	* 						first - total distance
	* 						second - number of crossed noed to get to "i" from source
	* 						this distance is a combination between the matching numbers between 2 faces and the 
	* deviation of the cross field (across the centers of all visited faces)
	* \param[out] previous     the id of the previous visited node
	* \param[out] found        the first singularity node we arrive into
	* \param[in] newNode2NewNodeNeighbours   the neighbours of a new node:
	*                          -  neigh(originalNode) - its original adjacent nodes + the new nodes corresponding 
	*						to the centers of its original adjacent faces
	* 					  -  neigh(newNode) - the original nodes of the face (whose center is the newNode)  + the new 
	* 						nodes corresponding to the centers of the original adjacent(be egde) faces 
	* \param[in] face2FaceTransport     the matching number between a face and its neighbours (same order as face2FaceNeighbours) 
	* \param[in] face2FaceDeviation     the deviation between a face and its neighbours (same order as face2FaceNeighbours) (calculated as the deviation of the cross field in between triangle centers)
	* \param[in] isBdryNode       boolean value indicating if the node is on the boundary  
	* \param[in] prevDir          previous direction
	* \param[in] prevDirCross     previous direction and previous cross (dictated by the previous node)
	* \param[in] targetDirCross   target (opposite) direction and cross (dictated by the node in which we want to arrive)
	* \param[in] forbiddenNodes boolean value indicating if we are allowed to pass through a node (we shouldn't be allowed to pass through other slot faces or singularity faces, except the source and target)
    	*/
	int getShortestPathBtwNewNodes(unsigned int&                                             numberOfNewNodes, 
	                               gmds::TCellID&                                            source,
	                               vector<unsigned int>&                                     targets,                                            
	                               vector<pair<double, unsigned int>>&                       min_distance,
	                               vector<int>&                                              previous, 
	                               int&                                                      found,
	                               vector<vector<gmds::TCellID>>&                            newNode2NewNodeNeighbours,
	                               vector<vector<double>>&                                   face2FaceTransport, 
	                               vector<vector<unsigned int>>&                             face2FaceDeviation,
	                               pair<gmds::math::Vector3d, gmds::math::Vector3d>&         prevDirCross,
	                               vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>>& targetDirCross,
	                               vector<bool>&                                             forbiddenFaces,
	                               bool&                                                     targetBdry,
	                               gmds::math::Point&                                        startPoint,
	                               gmds::math::Point&                                        endPoint,
	                               double&                                                   maxDist);
	
	/*----------------------------------------------------------------------------*/
	/** \brief computes the info in between 2 faces
	*/
	void computeFace2FaceInfo();
	/*----------------------------------------------------------------------------------------------------*/
	/* write a vtk file containing the singularity point locations and the slots*/
	void writeSingularityPointsAndSlots();
	/*----------------------------------------------------------------------------------------------------*/
	/*this strategy assumes that here we will impose the singularity triangles by providing their ID: Face singularTri
	   = m_mesh->get<Face>(triangle id); this has been implemented with the main scope of testing (ex for a cane-shaped
	   mesh, where the singular triangle is either not detected or detected in a different location that needed for
	   testing; note: the detection is highly dependent on the computation of the frame field; this strategy could also
	   be used if the user would have available a visualization platform and could simply click on a triangle -> it's ID
	   could be parsed directly)*/
	void initTestPrescribedSing();
	/*----------------------------------------------------------------------------------------------------*/
	void visualizeCrossVectors();
	/*----------------------------------------------------------------------------------------------------*/
	void deleteArtificialNodes(vector<CurveSingularityPoint *> &artificialSingPointsCreated);

	friend class LineRelaxator;
 protected:
    
	/** Mesh we start from */
	gmds::Mesh* m_mesh;
	/* Cross field we start from*/
	gmds::Variable<gmds::math::Cross2D>* m_field;
	/** DEBUG variable which stores the index of every single triangle of m_mesh */
	gmds::Variable<int> *m_index;

	Tools m_tool;
	/**directory where debug files will be written*/
	std::string m_output_directory_name;
    
	/** flag that indicates if geometric singularity points must be built */
	bool m_build_geometric_singularities;
	/** mark for nodes classified on geometric points */
	int m_mark_nodes_on_point;
	/** mark for nodes classified on geometric curves */
	int m_mark_nodes_on_curve;
	/** mark for edges classified on geometric curves */
	int m_mark_edges_on_curve;
    
	/**the obtained singularity graph*/
	SingularityGraph m_graph;
    
	/** radius of the bounding box containing m_mesh. */
	double m_mesh_radius;
	/** confusing distance used for connecting graph lines and points in the original strategy*/
	double m_confusing_distance;
    
	/*mean edge length of the original mesh*/
	double m_mean_edge_length;    
    
	double m_temp_epsilon; 
    
	unsigned int m_original_faces_number; 
    
	unsigned int m_original_nodes_number;
        
	bool m_withGlobalComments = false;
    
	/* \param m_ATolerance the tolerance used to connect singularity lines and
	* points. It is a % of the bounding box diagonal. It must be set in [0.01,0.1] */
	double m_ATolerance = 0.01; //default 0.01
    
	/** list of faces containing a 3-valent singularity point */
	std::list<gmds::Face> m_singularities_3;
	/** list of faces containing a 5-valent singularity point */
	std::list<gmds::Face> m_singularities_5;
    
	/** mark for faces containing a singularity point of the cross field*/
	int m_mark_faces_with_sing_point;
	/** mark for faces traversed by a singularity line of the cross field*/
	int m_mark_faces_with_sing_line;
    
	/** map used by the confusing ball. For faces, edges and nodes that are not
     located into a singularity, maps return value 0. */
	std::map<gmds::TCellID, SingularityPoint*> m_faces_to_singularity_on_surf;
	std::map<gmds::TCellID, SingularityPoint*> m_edges_to_singularity_on_surf;
	std::map<gmds::TCellID, SingularityPoint*> m_nodes_to_singularity_on_surf;   

	vector<gmds::math::Vector3d> m_face_normals;	
	vector<gmds::math::Point> m_triangle_centers; /* the points associated to the centers of all the triangles */
	vector<gmds::math::Cross2D> m_triangle_centers_cross; /*  the crosses associated to the centers of all the triangles */
};

// Manipulate Singularity Point locations from a grapgh in order to optimize the ratio of linked edges
class LineRelaxator
{
 public:
	LineRelaxator(SingularityGraph *graph, const double meanEdgeLength);

	void run();
	double getWorstAngle();
	double getWorstEdgeRatio();

 private:

	void computeFixedGroupAndFixedLine();
	void computeMeanEdgeLengthByGroup();
	void applySpringForceOnPoints();
	void movePoints();
	void MoveIfNewLocationIsWorthTheAnglePenalty(SingularityPoint *singPoint, const gmds::math::Point newLocation);

	SingularityGraph *m_graph;

	double m_step;
	std::vector<double> m_lineSpringLenght;
	std::vector<double> m_minLengthByGroup;
	std::vector<double> m_maxLengthByGroup;
	std::vector<SingularityLine *> m_minLineByGroup;
	std::vector<SingularityLine *> m_maxLineByGroup;

	std::vector<size_t> m_numberOfLinePerGroup;
	std::vector<bool> m_groupIsFixed;
	std::map<size_t, bool> m_lineIsfixed;
	std::map<size_t, gmds::math::Vector3d> m_directions;

	std::map < size_t, std::vector<double>> m_angleAtSingularity;
};
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYGRAPHBUILDER_H_ */
/*----------------------------------------------------------------------------*/
