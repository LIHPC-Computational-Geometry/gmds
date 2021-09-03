/*----------------------------------------------------------------------------*/
/*
 * ToolKit.h
 *
 *  Created on: sept. 18, 2015
 *      Author: Franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef TOOLS_H_
#define TOOLS_H_
/*----------------------------------------------------------------------------*/
// STL Headers
/*----------------------------------------------------------------------------*/
// gmds Headers
#include <gmds/ig/Mesh.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Cross2D.h>
#include <gmds/math/Chart.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Triangle.h>
#include <gmds/math/Segment.h>
#include <gmds/math/AxisAngleRotation.h>
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
/*----------------------------------------------------------------------------*/
/* \class This class gathers elementary algorithms used for different purposes
 *        by the main algorithms provided in FRAME.
 */
class LIB_GMDS_SINGGRAPHBUILD_API Tools
{
public:
    
    /**\brief Removes boundary slivers for a mesh with R|N|R2N|N2R model
     * \param[in,out] AMesh the mesh to work on
     */
    static int removeBoundarySlivers(gmds::Mesh* AMesh);
    
    struct PointVolumetricData{
        gmds::math::Point pnt;      // Point we are located at
        gmds::math::Vector3d dir;   // Direction we propage to
        gmds::Region tet;           // tetra where pnt is located in
        
        PointVolumetricData(const gmds::math::Point& AP,
                            const gmds::math::Vector3d& AV,
                            const gmds::Region& AR)
        :pnt(AP),dir(AV),tet(AR){;}
    };
    struct PointSurfacicData{
        gmds::math::Point pnt;    // Point we are located at
        gmds::math::Vector3d dir; // Direction we propage to
        gmds::Face tri;           // triangle where pnt is located in
        
        PointSurfacicData(const gmds::math::Point& AP,
                            const gmds::math::Vector3d& AV,
                            const gmds::Face& AF)
        :pnt(AP),dir(AV),tri(AF){;}
    };
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param AMesh the mesh where we work on
     * \param AField the cross field associated to AMesh
     */
    Tools(gmds::Mesh* AMesh,
          gmds::Variable<gmds::math::Cross2D>* AField,
          gmds::Variable<gmds::math::AxisAngleRotation>* ARotField=0);

    /*------------------------------------------------------------------------*/
    /** \brief  Compute where the frame field defined by m_axis_angle will
     *          transport \p AData.pnt considering data in \p AData and the max
     *          distance is \p AMaxDist.
     *
     * \details This algorithm uses Heun's scheme + a fuzzy approch where we 
     *          always go in/out of simplex along faces
     *
     * \param[in] AData     point data in the flow
     * \param[in] AMaxDist  Maximum distance allowed for the displacement
     *
     * \param[out] APnt the reached location.
     * 
     * \return true if APnt has been computed, false if we reach a tet 
     *              being FF singular
     */
    bool followFlow(const PointVolumetricData& AData,
                    const double               AMaxDist,
                    gmds::math::Point&         APnt);
    
    bool followFlow(const PointSurfacicData& AData,
                    const double             AMaxDist,
                    const int                AMarkEdgeOnCurve,
                    gmds::math::Point&       APnt);
    
    bool isFFSingular(const gmds::Region& AR);
    bool isFFSingular(const gmds::Face& AF);
    gmds::math::Chart::Mapping getRij(const gmds::TCellID AFrom,
                                      const gmds::TCellID ATo) const;
    
    /*------------------------------------------------------------------------*/
    /** \brief Gives the out point when we go through a triangle following the
     *         cross field
     *
     * \param AFace       the face we work on
     * \param AInPnt      the geometric point we start from
     * \param AInVec      the geometric direction to follow
     * \param AInCellDim  the dimension of the cell start_pnt is located
     * \param AInCellID   the id of the cell start_pnt is located on
     * \param AOutPnt     the geometric point we go out
     * \param AOutVec     the geometric direction to follow after
     * \param AOutCellDim the dimension of the cell we go out
     * \param AOutCellID  the id of the cell we go out
     */
    void	traverseTriangle(const gmds::Face&         AFace,
                             const gmds::math::Point&  AInPnt,
                             const gmds::math::Vector3d& AInVec,
                             const int                 AInCellDIm,
                             const gmds::TCellID       AInCellID,
                             gmds::math::Point&        AOutPnt,
                             gmds::math::Vector3d&       AOutVec,
                             int&                      AOutCellDIm,
                             gmds::TCellID&            AOutCellID,
					    double&          streamlineDeviation);
    
    /*------------------------------------------------------------------------*/
    /** \brief Gives the out point when we go through a triangle following the
     *         cross field and starting from a node
     *
     * \param AFace       the face we work on
     * \param ANode       the node we come from
     * \param AInPnt      the geometric point we start from
     * \param AInVec      the geometric direction to follow
     * \param AOutPnt     the geometric point we go out
     * \param AOutVec     the geometric direction to follow after
     * \param AOutCellDim the dimension of the cell we go out
     * \param AOutCellID  the id of the cell we go out
     */
    void	traverseTriangle(const gmds::Face&         AFace,
                             const gmds::Node&         ANode,
                             const gmds::math::Point&  AInPnt,
                             const gmds::math::Vector3d& AInVec,
                             gmds::math::Point&        AOutPnt,
                             gmds::math::Vector3d&       AOutVec,
                             int&                      AOutCellDIm,
                             gmds::TCellID&            AOutCellID,
					    double&           streamlineDeviation);
    
    /*------------------------------------------------------------------------*/
    /** \brief Gives the out point when we go through a triangle following the
     *         cross field and starting from an edge
     *
     * \param AFace       the face we work on
     * \param AEdge       the edge we come from
     * \param AInPnt      the geometric point we start from
     * \param AInVec      the geometric direction to follow
     * \param AOutPnt     the geometric point we go out
     * \param AOutVec     the geometric direction to follow after
     * \param AOutCellDim the dimension of the cell we go out
     * \param AOutCellID  the id of the cell we go out
     */
    void	traverseTriangle(const gmds::Face&         AFace,
                             const gmds::Edge&         AEdge,
                             const gmds::math::Point&  AInPnt,
                             const gmds::math::Vector3d& AInVec,
                             gmds::math::Point&        AOutPnt,
                             gmds::math::Vector3d&       AOutVec,
                             int&                      AOutCellDIm,
                             gmds::TCellID&            AOutCellID,
					    double&          streamlineDeviation);
    
    /*------------------------------------------------------------------------*/
    /** \brief Performs the Heun's algorithm into a triangle to compute the
     *                  out point when we go through this triangle and we arrive
     *                  from an edge. Note that the algorithm is numerically
     *                  corrected to avoid to go back in the face we come from.
     *
     * \param AInEdge   the edge we come from
     * \param AInPnt    the geometric point we start from
     * \param AInVec    the geometric direction to follow
     * \param AOppNode  the node opposite to AInEdge
     * \param AInNode1  the first node of AInEdge
     * \param AInNode2  the second node of AInEdge
     * \param AOutEdge1 another edge of the face we work on
     * \param AOutEdge2 a second another edge of the face we work on
     * \param AOutPnt   the geometric point we go out
     * \param AOutVec   the geometric direction to follow after
     *
     * \return an index indicating which is the cell we go out through:
     *         - 1 for AOppNode
     *         - 2 for AInNode1
     *         - 3 for AInNode2
     *         - 4 for AOutEdge1
     *         - 5 for AOutEdge2
     *
     */
    int heunsComputation(//const gmds::Edge&         AInEdge,
                         const gmds::math::Point&  AInPnt,
                         const gmds::math::Vector3d& AInVec,
                         const gmds::Node&         AOppNode,
                         const gmds::Node&         AInNode1,
                         const gmds::Node&         AInNode2,
                         const gmds::Edge&         AOutEdge1,
                         const gmds::Edge&         AOutEdge2,
                         gmds::math::Point&        AOutPnt,
                         gmds::math::Vector3d&       AOutVec,
					double&                   deviation);
    
    /*------------------------------------------------------------------------*/
    /** \brief Performs the Heun's algorithm into a triangle to compute the
     *                  out point when we go through this triangle and we arrive
     *                  from a node. Note that the algorithm is numerically
     *                  corrected to avoid to go out from the triangle we get
     *                  in.
     *
     * \param AInNode   the node we come from
     * \param AInPnt    the geometric point we start from
     * \param AInVec    the geometric direction to follow
     * \param AOppNode1 the first node of AOppEdge
     * \param AOppNode2 the second node of AOppEdge
     * \param AOppEdge  the edge opposite to AInNode
     * \param AOutPnt   the geometric point we go out
     * \param AOutVec   the geometric direction to follow after
     *
     * \return an index indicating which is the cell we go out through:
     *         - 0 if it does not intersect
     *         - 1 for AOppNode1
     *         - 2 for AOppNode2
     *         - 3 for AOppEdge
     *
     */
    int heunsComputation(//const gmds::Node&         AInNode,
                         const gmds::math::Point&  AInPnt,
                         const gmds::math::Vector3d& AInVec,
                         const gmds::Node&         AOppNode1,
                         const gmds::Node&         AOppNode2,
                         const gmds::Edge&         AOppEdge,
                         gmds::math::Point&        AOutPnt,		
                         gmds::math::Vector3d&       AOutVec,
				     double&                   deviation);
    
    
     /*------------------------------------------------------------------------*/
    /** \brief Performs the Runge-Kutta algorithm into a triangle to compute the
     *                  out point when we go through this triangle and we arrive
     *                  from an edge. Note that the algorithm is numerically
     *                  corrected to avoid to go back in the face we come from.
     *
     * \param AINPnt    the geometric point we start from
     * \param AINVec    the geometric direction to follow
     * \param[out] AOUTPnt   the geometric point we go out (in case we arrive at boundary -> the obtained bdry point)
     * \param[out] AOUTVec   the geometric direction to follow after
     * \param deviation measure of the deviation from the field
     * \param stepSize  the step size for the algorithm (∆t)
	* \param[out] AEndOnBdry  boolean value indicating if we have arrived to the boundary
     *\param[out] AToCellDim the dimension of the last visited element (node, edge, face)
	* \param[out] AToCellID the id of the last visited element (node, edge, face)
     */
    void RK4Computation(const gmds::math::Point&  AINPnt,
               const gmds::math::Vector3d& AINVec,             
               gmds::math::Point&        AOUTPnt,
               gmds::math::Vector3d&       AOUTVec,
               double&              deviation,
               double&              stepSize,	
			bool& AEndOnBdry,
			std::vector<gmds::TCellID>&     ATriangles,
			int&                            AToCellDim,
               gmds::TCellID&                  AToCellID);

    /*------------------------------------------------------------------------*/
    /** \brief Let a ray r=[AInPnt, AInVec) interesecting an edge AEdge, this
     *         method returns the point of intersection P between r and AEdge,
     *         and the output vector at P, that respects the underlying cross
     *         field
     *
     * \param AEdge     the edge we intersect
     * \param AInPnt    the geometric point we start from
     * \param AInVec    the geometric direction to follow
     * \param AOutPnt   the geometric point we go out
     * \param AOutVec   the geometric direction to follow after
     * \param deviation   deviation of the resulted line wrt the closest component of the cross at the intersection point
	* 
     * \return true if the ray intersects the edge, false otherwise
     */
    bool computeOutVectorFromRayAndEdge(const gmds::Edge&         AEdge,
                                        const gmds::math::Point&  AInPnt,
                                        const gmds::math::Vector3d& AInVec,
                                        gmds::math::Point&        AOutPnt,
                                        gmds::math::Vector3d&       AOutVec,
								double&             deviation);
    
    /*------------------------------------------------------------------------*/
    /** \brief Arriving at node ANode with direction AInVec, it returns the
     *         best fit direction prescribed by the cross field at ANode.
     *
     * \param ANode   the node we consider
     * \param AInVec    the geometric direction we arrive
     * \param AOutVec   the geometric direction to follow after
     *
     * \return true if the ray intersects the edge, false otherwise
     */
    void computeOutVectorAtPoint(const gmds::Node&         ANode,
                                 const gmds::math::Vector3d& AInVec,
                                 gmds::math::Vector3d&       AOutVec); 
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute the next cell we will go through. It is a face (dim=2),
     *         or an edge (dim=1).
     *
     * \param[IN]  AFromPnt     the geometric point we start from
     * \param[IN]  AFromVec     the geometric direction we go along
     * \param[IN]  AFromCellDim the dim. of the mesh cell we start from (0 or 1)
     * \param[IN]  AFromCellID  the id of the mesh cell we start from
     * \param[OUT] AToCellDim   the dim. of the mesh cell we start from (0 or 1)
     * \param[OUT] AToCellID    the id of the mesh cell we start from
     */
    void findNextCell(const gmds::math::Point&  AFromPnt,
                      const gmds::math::Vector3d& AFromVec,
                      const int AFromCellDim,
                      const gmds::TCellID AFromCellID,
                      int& AToCellDim,
                      gmds::TCellID& AToCellID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute the next cell we will got through. It is a face (dim=2),
     *         or an edge (dim=1).
     *
     * \param[IN]  AFromPnt   the geometric point we start from
     * \param[IN]  AFromVec   the geometric direction we go along
     * \param[IN]  AFromNode  the node we come from
     * \param[OUT] AToCellDim the dim. of the mesh cell we start from (0 or 1)
     * \param[OUT] AToCellID  the id of the mesh cell we start from
     */
    void findNextCell(const gmds::math::Point&  AFromPnt,
                      const gmds::math::Vector3d& AFromVec,
                      const gmds::Node& AFromNode,
                      int& AToCellDim,
                      gmds::TCellID& AToCellID);
    /*------------------------------------------------------------------------*/
    /** \brief Compute the next cell we will got through. It is a face (dim=2),
     *         or an edge (dim=1).
     *
     * \param[IN]  AFromPnt   the geometric point we start from
     * \param[IN]  AFromVec   the geometric direction we go along
     * \param[IN]  AFromEdge  the edge we come from
     * \param[OUT] AToCellDim the dim. of the mesh cell we start from (0 or 1)
     * \param[OUT] AToCellID  the id of the mesh cell we start from
     */
    void findNextCell(const gmds::math::Point&  AFromPnt,
                      const gmds::math::Vector3d& AFromVec,
                      const gmds::Edge& AFromEdge,
                      int& AToCellDim,
                      gmds::TCellID& AToCellID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Indicate if the ray starting from AFromNode and following AVec
     *         is aligned with AEdge.
     *
     * \param AVec      a direction modelized by a vector
     * \param AFromNode the node we start from
     * \param AEdge     an edge incident to AFromNode (otherwise returns false)
     *
     * \return a boolean indicating if the AVec is aligned with AEdge
     */
    bool isAlong(const gmds::math::Vector3d& AVec,
                 const gmds::Node& AFromNode,
                 gmds::Edge& AEdge);
    
    /*------------------------------------------------------------------------*/
    /** \brief Indicate if the ray starting from APnt and following AVec
     *         is going into the face AFace.
     *
     * \param APnt      the starting point
     * \param AVec      a direction modelized by a vector
     * \param AFromNode the node of AFace, that is located at APnt
     * \param AFace     the face we want to check
     *
     * \return a boolean
     */
    bool isGoingInto(const gmds::math::Point& APnt,
                     const gmds::math::Vector3d& AVec,
                     const gmds::Node& AFromNode,
                     const gmds::Face& AFace);
    /*------------------------------------------------------------------------*/
    /** \brief Indicate if the ray starting from APnt and following AVec
     *         is going into the face AFace.
     *
     * \param APnt      the starting point
     * \param AVec      a direction modelized by a vector
     * \param AFromEdge the edge of AFace, that contains APnt
     * \param AFace     the face we want to check
     *
     * \return a boolean
     */
    bool isGoingInto(const gmds::math::Point& APnt,
                     const gmds::math::Vector3d& AVec,
                     const gmds::Edge& AFromEdge,
                     const gmds::Face& AFace);
    
    gmds::math::Chart computeChartIn(const gmds::math::Point& APnt,
                                     const gmds::Face& AFace);
    
    void computeFuzzyHeuns(const gmds::math::Point&                 AFromPnt,
                           const gmds::math::Vector3d&              AFromDir,
                           const std::vector<gmds::Face>&           AFaces,
                           const std::vector<gmds::math::Triangle>& ATri,
                           gmds::math::Point&                       AToPnt,
                           gmds::math::Vector3d&                    AToDir,
                           int&                                     AToFaceId);
    void computeFuzzyHeuns(const gmds::math::Point&                 AFromPnt,
                           const gmds::math::Vector3d&              ADirPnt,
                           const std::vector<gmds::Edge>&           AFaces,
                           const std::vector<gmds::math::Segment>&  ATri,
                           gmds::math::Point&                       AToPnt,
                           gmds::math::Vector3d&                    AToDir,
                           int&                                     AToFaceId);
    
    /*------------------------------------------------------------------------*/
    /** \brief Indicates if a point \p APnt projected onto the plane defining
     *         \p ATri belong to \p ATri or not. Boolean \p AOnEdge0, \p AOnEdge1
     *         and \p \p AOnEdge2 indicates if the points lies on the edge o
     *         opposite to node 0, 1 and 3 respectively
     *
     * \param[in]  APnt a point
     * \param[in]  ATri a triangle
     * \param[out] AOnEdge0 true if APnt lies on the edge opposite to the node
     *                      0 of \p ATri
     * \param[out] AOnEdge1 true if APnt lies on the edge opposite to the node
     *                      1 of \p ATri
     * \param[out] AOnEdge2 true if APnt lies on the edge opposite to the node
     *                      2 of \p ATri
     *
     * \return true if \p APnt is in \p ATri, false otherwise
     */
    
    bool isPntInTri(const gmds::math::Point& APnt,
              const gmds::Face& ATri,
              bool& AOnEdge0, bool& AOnEdge1, bool& AOnEdge2, double& lambda1,
			  double& lambda2);
    
    
      /*------------------------------------------------------------------------*/
    /** \brief computes the next RK4 vector at a \p point_1 inside face \p AFace 
     * this function assumes we have already checked and \p point_1 is in triangle \p AFace
     *
     * \param[in]  point_1 a point
     * \param[in]  AFace a triangle
     * \param[in]  lambdas weights for barycentric coordinates of point \p point_1 inside triangle \p AFace
     * \param[in]  v_in the initial vector (whoose direction should be followed as close as possible by v_next)
     * \param[out] v_next the next RK4 vector
     */
    
    void getNextVectorRK4(//const gmds::math::Point& point_1,
                        const gmds::Face& AFace,
                        vector<double>& lambdas, 
                        const gmds::math::Vector3d& v_in,
                        gmds::math::Vector3d& v_next);
    
    
        /*------------------------------------------------------------------------*/
    /** \brief computes the next RK4 vector at a \p point_1 knowing that we have started from face \p AFace 
     * in this function we only know in which triangle(\p AFace) our initial point was located; therefore we have to find in which triangle \p point_1 is located
     *
	* \param[in]  AINPnt the original departing point
     * \param[in]  point_x the point at the current iteration
     * \param[in]  v_in the initial vector (whoose direction should be followed as close as possible by v_next)
     * \param[out] v_next the next RK4 vector
	* \param[in]  previousVisitedFaceq the triangles that have been visited previously (or the singularity triangle for the first computation)
	* \param[out] AEndOnBdry boolean value indicating if we have reached the boundary
	* \param[out] BdryPnt in case we have reached the boundary (AEndOnBdry == true), this is the final boundary point
	* \param[out] AToCellDim the dimension of the last visited element (node, edge, face)
	* \param[out] AToCellID the id of the last visited element (node, edge, face)
	* \param[out] isFinal boolean value indicating if we are athe last iteration of the RK4 step (in this case, AToCellDim and AToCellID must be modified to point to the actual element on whoch our last point lands)
     */
    
    void findTriangleAndNextVectorRK4(const gmds::math::Point& AINPnt,
							    gmds::math::Point& point_x,              
                        			   const gmds::math::Vector3d& v_in,
                       			   gmds::math::Vector3d& v_next,
							   bool& AEndOnBdry,
							   std::vector<gmds::TCellID>&     ATriangles,
							   int&                            AToCellDim,
               				   gmds::TCellID&                  AToCellID,
							   bool& isFinal,
							   vector<gmds::TCellID>& previousVisitedFaces);

private:
    
    /** Mesh we start from */
    gmds::Mesh* m_mesh;
    /* Cross field we start from*/
    gmds::Variable<gmds::math::Cross2D>* m_field;
    
    gmds::Variable<gmds::math::AxisAngleRotation>* m_rot_field;
    
};
/*----------------------------------------------------------------------------*/
#endif /* TOOL_H_ */
/*----------------------------------------------------------------------------*/
