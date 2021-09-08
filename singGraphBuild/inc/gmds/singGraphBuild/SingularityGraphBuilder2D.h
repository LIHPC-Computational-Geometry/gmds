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
	void execute(const Strategy AStrategy, unsigned int& number_of_control_points);
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
	void computeSingPointInfo(gmds::Face&        AFace,
	                          gmds::math::Point& APosSing);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief  For a a face AFace, which is singular, this method creates the
	* the singularity point and the slots (cont - used for visualization purposes)
	*/
	void createSingPointAndSlots(gmds::Face& AFace, unsigned int& cont);
    
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
	void buildGeometricSlots(vector<bool>& singOrGeomFaces);
	
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
	/** \brief Creation of singularity points
	*/
	//  void createSingularityPoints();
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creation of singularity lines - original strategy
	*/
	void createSingularityLines();
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creation of singularity lines starting simultaneously from all singularities 
	* (using Heun's computation)
	*/
	void createSingularityLinesSimultaneousStart();
	/*----------------------------------------------------------------------------------------------------*/
    
	/** \brief Creation of singularity lines starting simultaneously from all singularities 
	* (using RK4 computation)
	*/
	void createSingularityLinesSimultaneousStartRK4();
	/*----------------------------------------------------------------------------------------------------*/
    
	/** \brief Creation of singularity lines using the shortest paths method
	*/
	void createSingularityLinesShortestPaths(vector<vector<gmds::TCellID>>& newNode2NewNodeNeighbours,
	                                         vector<vector<double>>&        face2FaceTransport, 
	                                         vector<vector<unsigned int>>&  face2FaceDeviation,
	                                         vector<bool>&                  singOrGeomFaces);
	/*----------------------------------------------------------------------------------------------------*/
    
	/** \brief Remove a singularity line
	* \param ALine  the line that must be removed from the graph
	*/
	void removeSingularityLine(SingularityLine* ALine);
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Compute the singularity line starting from point ASingPnt in the
	*         direction of ASingSlot
	*
	* \param ASingPnt  the singularity point we start from
	* \param ASingSlot the associated slot we really consider
	* \param ARemovedSlot a slot that has been removed (if any)
	*/
	void computeSingularityLine(SingularityPoint*        ASingPnt,
	                            SingularityPoint::Slot*  ASingSlot, 
	                            unsigned int&            cont, 
	                            SingularityPoint::Slot*& ARemovedSlot);
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Compute the singularity line starting from point \param AFromPnt in the
	*         direction of \param AFromSlot
	*
	* \param[in] AFromPnt      the singularity point we start from
	* \param[in] AFromSlot     the associated slot we consider
	* \param[out] AToSingPnt    the singularity point we arrive at, if it exists
	* \param[out] AToSlot       the associated slot
	* \param[out] AToPnt        the location we arrive at the end
	* \param[out] AToDir        the last direction to arrive at the end
	* \param[out] APoints       a discretization of the stream line
	* \param[out] ATriangles    the list of traversed triangles
	* \param[out] AToCellDim    the dimension of the cell we finish on (0,1,2)
	* \param[out] AToCellID     the id of the cell we finish on
	* \param[out] streamlineDeviation     the deviation of the computed streamline
	* \param[out] AEndOnBnd     indicates if we finish on the boundary (true)
	* \param[out] AToSlotIsFree indicates if we finish onto a free slot (true)
	* \param[out] APntToCreate  indicates if we must create the end point (true)
	*/
	/*----------------------------------------------------------------------------------------------------*/
	void computeStreamLine(SingularityPoint*               AFromPnt,
	                       SingularityPoint::Slot*         AFromSlot,
	                       SingularityPoint*&              AToSingPnt,
	                       SingularityPoint::Slot*&        AToSlot,
	                       gmds::math::Point&              AToPnt,
	                       gmds::math::Vector3d&           AToDir,
	                       std::vector<gmds::math::Point>& APoints,
	                       std::vector<gmds::TCellID>&     ATriangles,
	                       int&                            AToCellDim,
	                       gmds::TCellID&                  AToCellID,
	                       double&                         streamlineDeviation,
	                       bool&                           AEndOnBnd,
	                       bool&                           AToSlotIsFree,
	                       bool&                           APntToCreate);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Improve a singularity line to avoid side effects due to the way
	*         singularity points are considered.
	*
	* \param ALine     the singularity line we work on
	* \param AFromPnt  the singularity point we start from
	* \param AFromSlot the associated slot we really consider
	* \param AToPnt    the singularity point we go to
	* \param AToSlot   the associated slot we really consider
	* \param to_dir    the direction of the singularity line
	* \param to_dir    the deviation of the newly computed singularity line
	* \param foundBackTrackPath   boolean value indicating wheather we have found a path
	*/
	void backtrackSingularityLine(SurfaceSingularityLine* ALine,
	                              SingularityPoint*       AFromPnt,
	                              SingularityPoint::Slot* AFromSlot,
	                              SingularityPoint*       AToPnt,
	                              SingularityPoint::Slot* AToSlot,
	                              gmds::math::Vector3d&   to_dir,
	                              bool&                   foundBackTrackPath);
    
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
	void growLineRK4(SingularityPoint*               AFromSingPnt,
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
	                 double&                         stepSize,
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
	void detectLineIntersections(const Strategy AStrategy);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief This operation creates intersections between 2 lines in the
	*         vicinity of a mesh face
	*/
	void createLineIntersection(const Strategy                  AStrategy,
	                            SurfaceSingularityLine*         ALine1,
	                            SurfaceSingularityLine*         ALine2,
	                            gmds::Face&                     AFace,
	                            std::vector<gmds::math::Point>& AAddedPoints);
	/*----------------------------------------------------------------------------------------------------*/
	void createLineIntersection(const Strategy AStrategy,
	                            std::vector<SurfaceSingularityLine*>& ALines,
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
	
/*----------------------------------------------------------------------------------------------------*/	

	/** \brief get the shortest path between a singular triangle \param[in] source 
	 * and all other singular triangles \param[in] targets by 'walking from face center to face center'
	\param[in] source    the singular triangle from which we start
	\param[in] targets   all the other singular triangles
	\param[in] targetPoints the end points of all slots (corresponding to \param[in] targets)
	\param[out] min_distance    min_distance[i] - minimum distance between source and face "i" - 
	* 						first - total distance
	* 						second - number of crossed faces to get to "i" from source.
	* 						this distance is a combination between the matching numbers between 2 faces and the 
	* deviation of the cross field (across the centers of all visited faces)
	\param[out] previous     the id of the previous visited face (previous[v] = u  mean that in our path the 
	shortest path detected from the \param source towards face "v" has the ending path={...u,v})
	\param[out] found        the first singularity face we arrive into
	\param[in] face2FaceMatch     the matching number between a face and its neighbours (same order as face2FaceNeighbours) 
	\param[in] face2FaceDeviation     the deviation between a face and its neighbours (same order as face2FaceNeighbours) (calculated as the deviation of the cross field in between triangle centers)
	\param[in] prevDirCrossGlobal     previous direction (first) and previous cross (second) - dictated by the previous triangle 
	\param[in] targetDirCross   the (opposite) slot direction for targets (in the same order as \param[in] targets)     
	\param[in] is_bdry_face          boolean value indicating whether a triangle is bounbdary (triangles having at least one node on the boundary are considered boundary triangles)
	\param[in] startPoint last added point (using RK4) for the begining slot; it is located inside face \param[in] source
	\param[out] finalEndPoint[i] detected point at the intersection between the boundary and the shortest path between source and the boundary
	\param[in] maxDist  the "infinity" value with which all weights are initialized; also used for penalizing paths that switch cross components
	\param[out] final_to_cell_dim last visited cell dimension - for boundary path
	\param[out] final_to_cell_id last visited cell id - for boundary path
	\param[out] finalLastVisTri last visited Face - for boundary path
	\param[in] contSource   value indicating the "order" of the source -> contSource = 5*i+j, where i represents the singularity point number and j represents the number of the slot (j<=5)
	\param[in] colNo   colNo = totalNumberOfSlots + 1; if we imagine the system as a matrix, the lines would be all the slots and the columns would be all the slots plus one column for boundary
	\param[in] totalNumberOfSlots   totalNumberOfSlots is actually the total number of possible slots: totalNumberOfSlots = 5 *|singularity points|
	\param[in] faceNo2Cont faceNo2Cont[triangle id] = 5*i+j, where i represents the singularity number and j represents the number of the slot (j<=5)
	\param[in] line_discretization the initial points (RK4) that belong to the \param[in] source initial line 
	\param[in] line_triangles the initial triangles that have been visited at the intial (RK4) step for \param[in] source  
	\param[out] pointPaths .first -> line_discretization for \param[in] source ; .second -> line_discretization for \param[in] targets[i]; this is needed only for visualization purposes 
	\param[out] finalPaths the list of shortest paths triangles (tri centers) from \param[in] source to \param[in] targets[i].
	\param[out] finalCenterPoints - the entire line_discretization associated to each shortest path line (line_discretization[source] + triangle_centers[finalPaths] + inverse (line_discretization[targets[i]]))
	\param[out] traversedSPTris - the entire set of visited triangles for each shortest path line
	NOTE: since we "walk" from triangle to neighbouring triangle by VERTEX, in between those two, we can have additional triangles that are detected using getTraversedTrisFaceNeighbyVertsExhaustive()
	\param[out] distances - weights for the detected shortest paths
	\param IllegalCross - pairs of mutually excluding shortest paths (that shouldn't co-exist in the final 
	graph)
	\param isTraversedNodeCompVector - stores for each node the list of shortest path "ids" that pass through it following the same cross component (first - horizontal, second - vertical)
	\param isTraversedFaceCompVector - stores for each face the list of shortest path "ids" that pass through it following the same cross component (first - horizontal, second - vertical)
	
	\param isTraversedFaceCompVector_SegmentPathCont - explained at \param getTraversedTrisFaceNeighbyVertsExhaustive
	\param isTraversedFaceCompVector_SegmentPathCode - explained at \param getTraversedTrisFaceNeighbyVertsExhaustive
	
	\param[in] contSourceToSingularity   for optimization; variable at location contMatrix( = contSource*colNo + i; if bdry => contSource*colNo+totalNumberOfSlots), 
	we store the singularity ids; at this location we have a sing line in between singularities i and j
			contSourceToSingularity[contMatrix].first = i;
			contSourceToSingularity[contMatrix].second = j;
	\param[in] singOrGeomFaces boolean value indicating whether a face is a singular face or a geometric one
	*/
	void getShortestPathBtwFacesOptimized(gmds::TCellID&                                                              source,
	                                      vector<unsigned int>&                                                       targets,
	                                      vector<gmds::math::Point>&                                                  targetPoints,
	                                      vector<pair<double, unsigned int>>&                                         min_distance, 
	                                      vector<int>&                                                                previous, 
	                                      int&                                                                        found,
	                                      vector<vector<double>>&                                                     face2FaceMatch, 
	                                      vector<vector<unsigned int>>&                                               face2FaceDeviation,
	                                      pair<gmds::math::Vector3d, gmds::math::Vector3d>&                           prevDirCrossGlobal,
	                                      vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>>&                   targetDirCross,
	                                      vector<bool>&                                                               is_bdry_face,
	                                      gmds::math::Point&                                                          startPoint,
	                                      vector<gmds::math::Point>&                                                  finalEndPoint,
	                                      double&                                                                     maxDist,
	                                      int&                                                                        final_to_cell_dim, 
	                                      gmds::TCellID&                                                              final_to_cell_id,
	                                      gmds::TCellID&                                                              finalLastVisTri,
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
	                                      vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedNodeCompVector,
	                                      vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedFaceCompVector,
	                                      vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedFaceCompVector_SegmentPathCode,
	                                      vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedFaceCompVector_SegmentPathCont,
	                                      vector<pair<unsigned int, unsigned int>>&                                   contSourceToSingularity,
	                                      vector<bool>&                                                               singOrGeomFaces);
		
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
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief retrace the shortest path between the \param[in] foundTarget and the source (from previously called getShortestPathBtwFaces)
	* \param[in] foundTarget     the target that we have reached previously (getShortestPathBtwFaces)
	* \param[in] previous        the vector storing for each position (corresponding to a face id) the previous visited face id
	* \param[out] path            the final path (from source to foundTarget)
	*/
  
	void retraceShortestPath(int&                  foundTarget,
	                         vector<int>&          previous, 
	                         vector<unsigned int>& path);
    
	void retraceShortestPath(gmds::TCellID&        foundTarget,
	                         vector<int>&          previous, 
	                         vector<unsigned int>& path);

	void retraceShortestPath(gmds::TCellID&        foundTarget,
	                         vector<int>&          previous, 
	                         vector<unsigned int>& path,
	                         unsigned int&         tempFacesNo);
	
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief having the shortest path between the target and the source as the sequence of triangles, detect the entire 
	* set of traversed triangles (since for the previous step we take into account neighboursByVert, the set of visited triangles 
	* obtained previously is not complete);
	* also, detect illegal crossing between potential singularity lines (if their corresponding cross component 
	* directions are alligned, they are not BOTH allowed to appear in the final graph)
	\param[in] line_discretization      the detected begining points for all slots (RK4)
	\param[in] line_triangles           the triangles that have been visited by line_discretization
	\param[in] finalPaths        			the vector containing path found (by retraceShortestPath) - the ids of the triangles in the shortest path
	\param[out] traversedSPTris      		the vector containing all the traversed triangles by shortest path line (from slot point to slot point)
	\param isTraversedNodeCompVector - stores for each node the list of shortest path "ids" that pass through it following the same cross component (first - horizontal, second - vertical)
	\param isTraversedFaceCompVector - stores for each face the list of shortest path "ids" that pass through it following the same cross component (first - horizontal, second - vertical)
	\param isTraversedFaceCompVector_SegmentPathCont - same order as in \param isTraversedFaceCompVector ; except that for each face it stores the actual segment that passes througn the face 
	"segment" stored as the position; the current segment will be segment(path[cont], path[cont+1]); we will store \param cont
	\param isTraversedFaceCompVector_SegmentPathCode - same order as in \param isTraversedFaceCompVector ; 
	since all paths are actually formed by source points + triangle centers (+target points), for each face it stores a code:
	 				0 - segment between 2 source points
	 				1 - segment between last source point and first triangle center
	 				2 - segment between 2 triangle centers(finalPaths[contSource][contTarget])
	 				3 - segment between the last triangle center and the first target point
	 				4 - segment between 2 target points
	 				5 - segment at boundary (last triangle center(finalPaths[contSource][contTarget]) and the finalEndPoint[contSource])
	 
	\param[in] contMatrix   value indicating the "order" of the path inside the optimization matrix -> contMatrix = contSource*colNo + totalNumberOfSlots;
	\param[in] colNo,  total number of columns of the system matrix -> totalNumberOfSlots+1 (1 corresponding to the boundary)
	\param[in] contSource   value indicating the "order" of the source -> contSource = 5*i+j, where i represents the singularity number and j represents the number of the slot (j<=5)
	\param[in] contTarget   value indicating the "order" of the target -> contSource = 5*i+j, where i represents the singularity number and j represents the number of the slot (j<=5)	 
	if boundary path, contTarget = totalNumberOfSlots = singPointNo*5;
	\param IllegalCross           		vector of pairs of pottential singularity lines that cannot co-exist in the final graph
	\param[in] contSourceToSingularity  unsigned int contMatrix = contSource*colNo + i;										
	\param[in] lastVisTri                    the last visited triangle - necessary only in the case of a bdry path	
	\param[in] finalEndPoint                    the last added points to the singularity lines - necessary only in the case of a bdry path
	*/
	
	void getTraversedTrisFaceNeighbyVertsExhaustive(vector<vector<vector<gmds::math::Point>>>&                  line_discretization,
	                                                vector<vector<vector<gmds::TCellID>>>&                      line_triangles,
	                                                vector<vector<vector<gmds::TCellID>>>&                      finalPaths,
	                                                bool&                                                       isTargetBdry,
	                                                vector<gmds::TCellID>&                                      traversedTriangles, 
	                                                vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedNodeCompVector,
	                                                vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector,
	                                                vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector_SegmentPathCode,
	                                                vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector_SegmentPathCont,
	                                                unsigned int&                                               contMatrix,
	                                                unsigned int&                                               colNo, //totalNumberOfSlots+1
	                                                unsigned int&                                               contSource,
	                                                unsigned int&                                               contTarget,
	                                                vector<pair<unsigned int, unsigned int>>&                   IllegalCross,
	                                                vector<pair<unsigned int, unsigned int>>&                   contSourceToSingularity,
	                                                gmds::TCellID&                                              lastVisTri,
	                                                vector<gmds::math::Point>&                                  finalEndPoint,
	                                                unsigned int&                                               i_orig,
	                                                unsigned int&                                               j_orig,
	                                                unsigned int&                                               t1_orig,
	                                                unsigned int&                                               t2_orig);   
	
	/*----------------------------------------------------------------------------------------------------*/
	
	/** \brief having the traversed triangles, detect the illegal crossing between potential singularity lines (if their corresponding cross component 
	* directions are alligned, they are not BOTH allowed to appear in the final graph)   
	\param[in] discrPtsSource        			the vector all the points of the potential singularity line (as detected by growLineRK4)
	\param[in] travTriSource                      the vector containing all the visited triangles of the potential singularity line (as detected by growLineRK4)
	\param[out] traversedTriangles      		the vector containing all the traversed triangles - this should technically be identical to travTriSource //WARNING TODO check!
	\param[out] isTraversedFace  			vector of boolean values indicating whether the triangles have been visisted by the temp singularity lines  (SP)
	\param[in] contMatrix         	  indicates the position of the current singularity line in the matrix 									
	\param[in] colNo                         the number of columns of our matrix (totalNumberOfSlots+1) -> 1 is for the bdry	
	\param[out] IllegalCross           		vector of pairs of pottential singularity lines that cannot co-exist in the final graph
	\param[in] contMatrixToSingularity           for each singularity line treated (treated in contMatrix order) store its corresponding singularity slot (the order 5*i+j, where i = sing number, j = slot number)
	\param[in] lastVisTri                    the last visited triangle - necessary only in the case of a bdry path
	\param[in] endPoint                    the last added point to the singularity line - necessary only in the case of a bdry path
	\param[in] illegalProblematicBdryPair  indicates the pair (in IllegalCross) of the current singularity line; technically we would want that the current singularity line is no longer paired with illegalProblematicBdryPair 
	\param[out] recomputeIllegalProblematicBdryPair  if even after the recomputation of the current singulaity line (using growLineRK4) this line still intersects in an illegal manner with illegalProblematicBdryPair, 
	\param[in] reAddIt                              if our problematic pair is stil found as illegal then this argument (reAddIt) indicates whether we should re-add it to the vector IllegalCross
	*/
	void getIllegalCrossByTraversedTris(vector<vector<vector<gmds::math::Point>>>&                  line_discretization,
	                                    vector<vector<vector<gmds::TCellID>>>&                      line_triangles,
	                                    vector<gmds::TCellID>&                                      traversedTriangles,
	                                    vector<vector<vector<gmds::TCellID>>>&                      finalPaths,
	                                    vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedNodeCompVector,
	                                    vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector,
	                                    vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector_SegmentPathCode,
	                                    vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector_SegmentPathCont,
	                                    unsigned int&                                               contMatrix,
	                                    unsigned int&                                               colNo, //totalNumberOfSlots+1
	                                    unsigned int&                                               contSource,
	                                    unsigned int&                                               contTarget,
	                                    vector<pair<unsigned int, unsigned int>>&                   IllegalCross,
	                                    vector<pair<unsigned int, unsigned int>>&                   contMatrixToSingularity,
	                                    gmds::TCellID&                                              lastVisTri,
	                                    vector<gmds::math::Point>&                                  finalEndPoint,
	                                    unsigned int&                                               illegalProblematicBdryPair,
	                                    bool&                                                       recomputeIllegalProblematicBdryPair,
	                                    bool&                                                       reAddIt,
	                                    unsigned int&                                               i_orig,
	                                    unsigned int&                                               j_orig);
   
	/*----------------------------------------------------------------------------*/
	/** \brief computes the info in between 2 faces
	* \param face2FaceTransport     the matching number between a face and its neighbours (same order as face2Face_neighbours) 
	* \param face2FaceDeviation     the deviation between a face and its neighbours (same order as face2Face_neighbours) (calculated as the deviation of the cross field in between triangle centers)
	*/
	void computeFace2FaceInfo(vector<vector<gmds::TCellID>>& newNode2NewNodeNeighbours,
	                          vector<vector<double>>&        face2FaceTransport, 
	                          vector<vector<unsigned int>>&  face2FaceDeviation);
	/*----------------------------------------------------------------------------------------------------*/
	/* write a vtk file containing the singularity point locations and the slots*/
	void writeSingularityPointsAndSlots();

	/* detect the list of faces surounding a vertex and storing it in vertNeighFaces*/
	void getVertNeigh();
    
	/*----------------------------------------------------------------------------------------------------*/

	
private:
    
	/** Mesh we start from */
	gmds::Mesh* m_mesh;
	/* Cross field we start from*/
	gmds::Variable<gmds::math::Cross2D>* m_field;
    
	/** DEBUG variable which stores the index of every single triangle of m_mesh */
	gmds::Variable<int>* m_index;
	/* the list of faces surounding a vertex */
	std::vector<std::vector<gmds::TCellID>> vertNeighFaces = std::vector<std::vector<gmds::TCellID>>(0, std::vector<gmds::TCellID>(0)); 
    
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
	/** technical container, which is used to store all the free slots during the
	*  algorithm */
	std::list<SingularityPoint::Slot*> m_free_slots;
    
	/** radius of the bounding box containing m_mesh. */
	double m_mesh_radius;
	/** confusing distance used for connecting graph lines and points in the original strategy*/
	double m_confusing_distance;
    
	/*mean edge length of the original mesh*/
	double mean_edge_length;    
    
	double temp_epsilon; 
    
	double original_faces_number; 
    
	double original_nodes_number;
        
	bool withGlobalComments = false;
	
	bool visualizeShortestPaths = true;
    
	/* \param ATolerance the tolerance used to connect singularity lines and
	* points. It is a % of the bounding box diagonal. It must be set in [0.01,0.1] */
	double ATolerance = 0.01; //default 0.01
    
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
    
	vector<gmds::math::Vector3d> face_normals;	
	vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> local_basis;/*local_basis calculated as considering the first edge as Ox axis and as the Oy axis as the cross product between Ox and the face normal*/
	vector<gmds::math::Vector3d> bdry_edge_normals;
	vector<gmds::math::Vector3d> bdry_node_normals;
	vector<vector<gmds::Face>> face2Face_neighbours; /*the neighbours of a face (considered as the adjacent faces to the 3 edges of the current face)*/
	vector<vector<gmds::Face>> face2Face_neighbours_by_verts; /*the neighbours of a face (considered as the adjacent faces to the 3 nodes of the current face)*/
	vector<bool> is_bdry_face; /*vector of boolean values indicating whether a triangle is situated along the boundary (true - even if only one of its nodes is on the boundary)*/
	vector<gmds::math::Point> triangle_centers; /* the points associated to the centers of all the triangles */
	vector<gmds::math::Cross2D> triangle_centers_cross; /*  the crosses associated to the centers of all the triangles */
    
	/*std::list<gmds::TCellID> m_2SingTetIDsAlongSurf;
	std::map<gmds::TCellID, int> m_region_singularity_type;
	std::map<gmds::TCellID, SurfaceSingularityPoint*> m_face_to_singularity_on_surf;
     
	//boolean marks for executing the algorithm
     
	int m_markNodesOnVert;
	int m_markNodesOnCurv;
	int m_markNodesOnSurf;
	int m_markAloneNodes;
     
	int m_markEdgesOnSurf;
	int m_markEdgesOnCurv;
     
	int m_markFacesOnSurf;
     
	int m_markClusterSingDone;
	int m_markBordVolSingForFaces;
	int m_markBordVolSingForEdges;
	int m_markBordVolSingForNodes;
     
	int m_markBordSurfSingForFaces;//faces containing sing. are marked with it
	int m_markBordSurfSingForEdges;//edges containing sing. are marked with it
	int m_markBordSurfSingForNodes;//nodes containing sing. are marked with it
	*/
};
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYGRAPHBUILDER_H_ */
/*----------------------------------------------------------------------------*/
