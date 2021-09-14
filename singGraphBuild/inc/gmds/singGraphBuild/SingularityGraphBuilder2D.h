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
#include <memory>
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
namespace gmds {

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
	SingularityGraphBuilder2D(Mesh* AMesh,
	                          Variable<math::Cross2D>* AField,
	                          const bool ABuildGeomSing = true);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Execution of the algorithm with the chosen AStrategy
	*/
	void execute();
	/*----------------------------------------------------------------------------------------------------*/
	/* \brief Returns a unique ptr of the mesh extracted fom the singularity graph 	*/
	std::unique_ptr<Mesh> getQuadMesh();
	/*----------------------------------------------------------------------------------------------------*/
	/* \brief Allows to optimize the final quad mesh : bring edge ratio to closer to 1.0 as well as angle closer to 90° */
	SingularityGraphBuilder2D& setQuadMeshSmoothingEnable(bool enableQuadMeshSmoothing);
	/*----------------------------------------------------------------------------------------------------*/
	/* \brief Optimizes the edges Ratio to be closed to 1 as well as angle close to 90° */
	SingularityGraphBuilder2D& setDebugFilesWritingEnable(bool enableDebugFilesWriting);
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
	* \param AMarkNodeForbiddenBdry mark for nodes on "forbidden curve" : 
				curves on which no singularity will be created
	*/
	void initMarks(const int AMarkNodePnt, const int AMarkNodeCrv,
	               const int AMarkEdgeCrv, const int AMarkNodeForbiddenBdry);

	/*----------------------------------------------------------------------------------------------------*/
    
protected:
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
	void computeSingPointInfo(const Face&        AFace,
	                          math::Point& APosSing);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief  For a a face AFace, which is singular, this method creates the
	* the singularity point and the slots (cont - used for visualization purposes)
	*/
	void createSingPointAndSlots(const Face& AFace);
    
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
	void writeTestMeshVerts(vector<Node>& ANodes, std::string& AFileName);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having as vertices \param APoints
	* and 'zero-area' triangles
	*/
	void writeTestPoints(vector<math::Point>& APoints, std::string& AFileName);
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having as triangles \param ATriangles
	*/
	void writeTestMeshTriangles(vector<Face>& ATriangles, std::string& AFileName);    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Write a new mesh named \param AFileName having as triangles a subset of the original mesh's triangles (\param ATrianglesIds), given as ids
	*/
	void writeTestMeshTrianglesIds(vector<TCellID>& ATrianglesIds, std::string& AFileName); 
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
	              math::Point&              AFromPnt,
	              math::Vector3d&           AToDir,
	              std::vector<math::Point>& APoints,
	              std::vector<TCellID>&     ATriangles,
	              int&                            AToCellDim,
	              TCellID&                  AToCellID,
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
	                 const math::Point&        AFromPnt,
	                 math::Vector3d&           AToDir,
	                 std::vector<math::Point>& APoints,
	                 std::vector<TCellID>&     ATriangles,
	                 int&                            AToCellDim,
	                 TCellID&                  AToCellID,
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
    
	void redefineOneConfusingBall(SingularityPoint* APnt, vector<TCellID>& modifiedFaces, double &previousRad, const double increaseRadiusScale = 1.5);
    
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
	void createGeometricSingularityPoint(const math::Point&    AInPnt,
	                                     const math::Vector3d& AInVec,
	                                     const int                   ACellDim,
	                                     const TCellID         ACellID,
	                                     SingularityPoint*&          APnt,
	                                     SingularityPoint::Slot*&    ASlot);
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Returs the singularity lines that go througt the face \param AFace
	*
	* \param AFace a mesh face
	* \return the singularity lines going through AFace
	*/
	std::vector<SurfaceSingularityLine*> getSingularityLinesIn(const Face& AFace);
    
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
	                            Face&                     AFace,
	                            std::vector<math::Point>& AAddedPoints);
	/*----------------------------------------------------------------------------------------------------*/
	void createLineIntersection(std::vector<SurfaceSingularityLine*>& ALines,
	                            Face&                           AFace,
	                            std::vector<math::Point>&       AAddedPoints);    
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Computes the barycenters for all triangles in the mesh, 
	* as well as the crosses associated to them
	*/
	void constructCenterTriangleCrosses();
    
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Computes the barycenter (\param centerTri) for triangle \param ATriangle, 
	*as well as the cross associated to it (\param centerTriCross) */
	void constructOneCenterTriangleCross(Face& ATriangle, 
	                                     math::Point &centerTri, 
	                                     math::Cross2D& centerTriCross);
    
	/*----------------------------------------------------------------------------------------------------*/
     
	/** \brief Constructs a new mesh (\param newLocalMesh) that has as Nodes all Nodes of 
	* \param trianglesToRemesh plus their centers and the resulted new triangles 
	* (each original triangle splits into 3 triangles)
	* \param newLocalMesh the newly constructed mesh    
	* \param newLocalMesh_id_to_mesh_id_node in this vector at the position corresponding the Node id of the \param newLocalMesh (except new introduced vertices) it is stored the id of the corresponding vertex in the original mesh 
	* \param local_cross_field_2D cross field of the \param newLocalMesh
	* \param trianglesToRemesh the ids of the triangles of the original mesh that we want to remesh
	*/
	void remeshTriangles(Mesh*                          newLocalMesh,  
	                     vector<TCellID>&               newLocalMesh_id_to_mesh_id_node, 
	                     Variable<math::Cross2D>* local_cross_field_2D, 
	                     vector<TCellID>&               trianglesToRemesh);
     
     
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
	void remeshTrianglesNewMesh(Mesh*                          newLocalMesh,   
	                            Variable<math::Cross2D>* newMesh_cross_field_2D, 
	                            vector<TCellID>&               trianglesToRemesh,
	                            vector<bool>&                        remeshedTriangles,
	                            vector<math::Point>&           newTriangleCenters, 
	                            vector<math::Cross2D>&         newTriangleCenterCrosses, 
	                            vector<math::Vector3d>&        newBdryEdgeNormals,
	                            vector<math::Vector3d>&        newBdryNodeNormals,
	                            vector<bool>&                        isCurveEdge,
	                            vector<bool>&                        isCurveNode,
	                            vector<TCellID>&               newLocalMesh_id_to_mesh_id_node,
	                            vector<TCellID>&               mesh_id_to_newLocalMesh_id_node,
	                            vector<TCellID>&               newLocalMesh_id_to_mesh_id_face,
	                            vector<TCellID>&               mesh_id_to_newLocalMesh_id_face);
    
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
     
	void connectSingularityLines(math::Vector3d            dir_slot_i,  
	                             math::Vector3d            dir_slot_j,
	                             math::Vector3d            connection_dir,
	                             math::Point               start_pnt1,// start_pnt = line_discretization[0]; //starting point
	                             math::Point               start_pnt2, 
	                             TCellID                   start_cell_id1,
	                             TCellID                   start_cell_id2,
	                             int                             start_cell_dim1, 
	                             int                             start_cell_dim2,
	                             std::vector<math::Point>& APoints,
	                             std::vector<TCellID>&     ATriangles,                 
	                             double&                         streamlineDeviation);
	
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
	/*----------------------------------------------------------------------------------------------------*/
	double computeMeanEdgeLength();

 protected:
    
	/** Mesh we start from */
	Mesh* m_mesh;
	/* Cross field we start from*/
	Variable<math::Cross2D>* m_field;

	Tools m_tool;
	/**directory where debug files will be written*/
	std::string m_output_directory_name;
    
	/** flag that indicates if geometric singularity points must be built */
	bool m_build_geometric_singularities;

	// @ { input marks provided by user
	/* mark for nodes classified on geometric points */
	int m_mark_nodes_on_point;
	/* mark for nodes classified on geometric curves */
	int m_mark_nodes_on_curve;
	/* mark for edges classified on geometric curves */
	int m_mark_edges_on_curve;
	 /* mark for forbidden boundary node */
	int m_mark_forbiddenBdryNode;
	// @ } 

	// @ {  internal marks
    /* mark on forbidden boundary edge*/
	int m_mark_forbiddenBdryEdge;
	/* mark for faces containing a singularity point of the cross field*/
	int m_mark_faces_with_sing_point;
	/* mark for faces traversed by a singularity line of the cross field*/
	int m_mark_faces_with_sing_line;
    // @ } 

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
        
	bool m_enableQuadMeshSmoothing = false;
	bool m_enableDebugFilesWriting = false;

	bool m_withGlobalComments = false;
    
	/* \param m_ATolerance the tolerance used to connect singularity lines and
	* points. It is a % of the bounding box diagonal. It must be set in [0.01,0.1] */
	double m_ATolerance = 0.01; //default 0.01
    
	/** list of faces containing a 3-valent singularity point */
	std::list<Face> m_singularities_3;
	/** list of faces containing a 5-valent singularity point */
	std::list<Face> m_singularities_5;
        
	/** map used by the confusing ball. For faces, edges and nodes that are not
     located into a singularity, maps return value 0. */
	std::map<TCellID, SingularityPoint*> m_faces_to_singularity_on_surf;
	std::map<TCellID, SingularityPoint*> m_edges_to_singularity_on_surf;
	std::map<TCellID, SingularityPoint*> m_nodes_to_singularity_on_surf;   

	vector<math::Vector3d> m_face_normals;	
	vector<math::Point> m_triangle_centers; /* the points associated to the centers of all the triangles */
	vector<math::Cross2D> m_triangle_centers_cross; /*  the crosses associated to the centers of all the triangles */
};


}
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYGRAPHBUILDER_H_ */
/*----------------------------------------------------------------------------*/
