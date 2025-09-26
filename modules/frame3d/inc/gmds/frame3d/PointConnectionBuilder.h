/*----------------------------------------------------------------------------*/
#ifndef GMDS_FRAME3D_GRID_POINT_CONNECTOR_H_
#define GMDS_FRAME3D_GRID_POINT_CONNECTOR_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <gmds/math/AxisAngleRotation.h>
#include <gmds/utils/OrientedGraph.h>
#include "GMDSFrame3d_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  PointConnectionBuilder
 *  \brief  This class provides an algorithm to connect a set of points that are
 *          defined on tetrahedral mesh where a frame field lives on. A frame
 *          is defined for each point to connect. The output is a set of "connection"
 *          information indicating for each frame field direction where a point
 *          is connected to.
 */
class GMDSFrame3d_API PointConnectionBuilder{

public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] AMesh   background tetrahedral mesh \p APnts was built
     *                    from
     * \param[in] APnts   the points we build the hex-dom mesh from
     * \param[in] ACharts the charts associated to each point in \p APnts
     * \param[in] ATypes  a flag indicating if the corresponding point was
     *                    extracted from a stable tet (0), a PGP sing tet (1)
     *                    or a FF sing tet(2)
     * \param[in] AClass  a flag indicating if the corresponding point is
     *                    classified onto a point (0), a curve (1), a surf(2)
     *                    or a volume (3)
     * \param[in] ACurv   a flag indicating the curve number it is classified
     *                    with 0 value if it is in the volume, on a surface or a
     *                    point
     * \param[in] ASurf   a flag indicating the surface number it is classified
     *                    with 0 value if it is in the volume, on a curve or a
     *                    point
     * \param[in] ANormal normal for each boundary node
     */
    PointConnectionBuilder(gmds::Mesh* AMesh,
                    const std::vector<gmds::math::Point>& APnts,
                    const std::vector<gmds::math::Chart>& ACharts,
                    const std::vector<gmds::Cell::Data>&  AData,
                    const std::vector<int>&               ATypes,
                    const std::vector<int>&               AClass,
                    const std::vector<int>&               ACurv,
                    const std::vector<int>&               ASurf,
                    const std::map<int, gmds::math::Vector3d>& ANormal);

    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for generating the mesh
     */
    void execute();
    /*------------------------------------------------------------------------*/
    /** \brief Returns edge info. Once performed the algorithm, we can acces to
     *         edges which are defined by the point index int the vector of
     *         point given when constructing *this
     *
     * @param[out] AEdges built edges.
     */
    void getEdges(std::vector<std::pair<int,int > >& AEdges);
    /*------------------------------------------------------------------------*/
    /** \brief Returns hexa info. Once performed the algorithm, we can acces to
     *         all the hexahedral cells that can be created
     *
     * @param[out] AHex built hexes
     */
    void getHexes(std::vector<std::vector<int> >& AHexes);

    /*------------------------------------------------------------------------*/
    /** \brief Access to the generated mesh
     *
     * \return the hexahedral cells we are able to generate
     */
    const gmds::Mesh& mesh() const { return m_hexes; }
    /*------------------------------------------------------------------------*/
    /** @brief set the debug info flag used during the algorithm execution.
     *
     * @param AWithDebug true to output debug info (log and files)
     * @param AOutputDir directory where debug files should be written in
     */
    void setDebugInfo(const bool& AWithDebug,
                      const std::string& AOutputDir=".");

    /*------------------------------------------------------------------------*/
    /** @brief In order to minimise some spatial research a spacing parameter
     *         is used in the code. Its setting is not obvious, so we can change
     *         it by this way for experiments
     *
     * @param ASpacing spacing value used in spatial search
     */
    void setSpacing(const double ASize=1.0){m_spacing=ASize;}
protected:

    enum EPointClassification{
        ON_VERTEX =0,
        ON_CURVE  =1,
        ON_SURFACE=2,
        IN_VOLUME =3
    };

    enum EPointType{
        REGULAR    =0,
        PARAM_SING =1,
        FRAME_SING =2
    };

    /*------------------------------------------------------------------------*/
    /** \struct HexCorner
     * \brief Local Structure to store hex corner data for each input point. A
     *        corner data is made of the index of the current point (p), the
     *        index of the 3 points helping to form the corner (adj). Each
     *        "edge" [p,adj[i]] is geometrically defined by vec[i].
     *
     *        The free field indicates if the corner is already used for a hex.
     *        The index field give the position of *this in the list of corners
     *        defined in point this->p
     *
     */
    struct HexCorner {
        int p;
        int index;
        bool free;
        int adj[3];
        gmds::math::Vector3d vec[3];
    };
    /*------------------------------------------------------------------------*/
    /** \struct OrientedEdge
     * \brief Local Structure to store an oriented edge made of the index of
     *        the first and second points and the chart data used in the first
     *        point to select this OrientedEdge. A chart data (-1,-1) means
     *        that the selected end point finally do not correspond to a chart
     *        alignment (will happen for boundary points for instance)
     */
    struct OrientedEdge {
        int first;
        int second;
        int axis;
        int dir;
        explicit OrientedEdge(const int AF=-1, const int AS=-1,
                              const int AAxis=-1, const int ADir=-1):
                first(AF),second(AS),axis(AAxis), dir(ADir){}

        bool operator==(const OrientedEdge& AE)const{
            return (first==AE.first) && (second==AE.second);
        }
        void setInvalid(){
            first=-1; second=-1;
        }
        bool isInvalid() const{
            return (first==-1 && second==-1);

        }

        bool isInverse(OrientedEdge& AE) const{
            return (AE.first==second && AE.second==first);
        }
        bool isEqual(OrientedEdge& AE) const{
            return (AE.first==first && AE.second==second);
        }

        /*--------------------------------------------------------------------*/
        /** \brief Weak equality (equal or inv) between  *this and \p AE
         *
         * \param[in] AE an oriented edge bo be compared with
         *
         * \return true if \p AE is the same or the inverse edge of *this,
         *         false otherwise.
         */

        bool operator==(OrientedEdge& AE) const{
            return isEqual(AE) || isInverse(AE);
        }

    };

    /*------------------------------------------------------------------------*/
    /** \brief Build for each point a list of close points to be compared with
     *         during edge creation.
     */
    void createDistanceFilter();


    /*------------------------------------------------------------------------*/
    /** \brief Associate each point to a mesh cell
     */
    void computeMeshAssociation();


    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges from each point to its neighbor points in a
     *         local point of view. \p AEdges[i] contains all the oriented edges
     *         coming from point i.
     *
     * \param[out] AEdges the list of oriented edges that will be created.
     */
    void buildOrientedEdges(std::vector<std::vector<OrientedEdge> >& AEdges);

    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges for the point \p APntId which is located
     *         into the volume. The input point is necessary REGULAR.
     *
     * \param[in]  APntID index of the point we work on
     * \param[in]  AVec   the [3][2] vectors of the frame in point APntID
     * \param[out] APnt   the [3][2] potential points location found
     * \param[out] AIndex the [3][2] potential point  index found
     * \param[out] AFound the [3][2] Indicate if a candidate is found or not
     */
    void buildOrientedEdgesInVolume(const int                  APntID,
                                    const gmds::math::Vector3d AVec[][2],
                                    gmds::math::Point          APnt[][2],
                                    int                        AIndex[][2],
                                    bool                       AFound[][2]);

    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges for the point \p APntId which is located
     *         onto a surface. The input point is necessary REGULAR.
     *
     *
     * \param[in]  APntID index of the point we work on
     * \param[in]  AVec   the [3][2] vectors of the frame in point APntID
     * \param[out] APnt   the [3][2] potential points location found
     * \param[out] AIndex the [3][2] potential point  index found
     * \param[out] AFound the [3][2] Indicate if a candidate is found or not
     */
    void buildOrientedEdgesOnSurface(const int                  APntID,
                                     const gmds::math::Vector3d AVec[][2],
                                     gmds::math::Point          APnt[][2],
                                     int                        AIndex[][2],
                                     bool                       AFound[][2]);

    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges for the point \p APntId which is located
     *         on a curve. The input point is necessary REGULAR.
     *
     * \param[in]  APntID index of the point we work on
     * \param[in]  AVec   the [3][2] vectors of the frame in point APntID
     * \param[out] APnt   the [3][2] potential points location found
     * \param[out] AIndex the [3][2] potential point  index found
     * \param[out] AFound the [3][2] Indicate if a candidate is found or not
     */
    void buildOrientedEdgesOnCurve(const int                  APntID,
                                   const gmds::math::Vector3d AVec[][2],
                                   gmds::math::Point          APnt[][2],
                                   int                        AIndex[][2],
                                   bool                       AFound[][2]);


    /*------------------------------------------------------------------------*/
    /** \brief Select among different point which is the bes aligned and/or
     *         at a good distance to build an oriented edge
     *
     * \param[in] AID   the candidate point we work
     * \param[in] ADot  the dot product of the ref pnt and the candidate pnt
     * \param[in] ADist the distance between the ref pnt and the candidate pnt
     *
     * \return    the index of the "best" candidate
     */

    int filterPointsForBuildingOrientedEdge(const std::vector<int>& AID,
                                            const std::vector<double>& ADot,
                                            const std::vector<double>& ADist);
    /*------------------------------------------------------------------------*/
    /** \brief Correct oriented edges to get a onsistent structure
     *
     * \param[in]  AInEdges  the list of oriented edges that we start from.
     * \param[out] AOutEdges the list of oriented edges that will be created.
     */
    void buildEdges(std::vector<std::vector<OrientedEdge> >& AInEdges,
                    std::vector<std::vector<OrientedEdge> >& AOutEdges);
    /*------------------------------------------------------------------------*/
    /** \brief ...
     */
    bool computeVolumePointFrom(const int APntIndex,
                                const gmds::math::Vector3d& AV,
                                gmds::math::Point& APnt);
    /*------------------------------------------------------------------------*/
    /** \brief ...
     */
    bool computeSurfacePointFrom(const int APntIndex,
                                 const gmds::math::Vector3d& AV,
                                 gmds::math::Point& APnt);

    /*------------------------------------------------------------------------*/
    /** \brief Check if an oriented edge going from \p AFrom to \p ATo exists
     *         in \p AEdgeSet.
     *
     * \param[in]  AFrom    A first point index
     * \param[in]  ATo      A second point index
     * \param[in]  AEdgeSet A set of oriented edges
     * \param[out] AOutEdge The edge from \p AFrom to \p ATo if it exists
     *
     * \return true if an edge going from \p AFrom to \p ATo exists in
     *         \p AEdgeSet, false otherwise.
     */
    bool isIn(const int AFrom, const int ATo,
              const std::vector<OrientedEdge>& AEdgeSet,
              OrientedEdge& AOutEdge) const;

    /*------------------------------------------------------------------------*/
    /** \brief Fill the m_hc_mapping attribute for each stable points. Such a
     *         point should be the corner of 1, 2, 4 or 8 hexahedral elts
     *         depending if it is classified on a geometric vertex, curve,
     *         surface or volume.
     */
    void buildHexCorners(std::vector<std::vector<OrientedEdge> >& AEdges);

    /*------------------------------------------------------------------------*/
    /** \brief Build hexahedral elements
     */
    void buildHexahedral();

    /*------------------------------------------------------------------------*/
    /** \brief Find among the free corners of \p AI and \p AJ an end point
     *         which is different of \p AFrom but both in a free corner of
     *         \p AI and \p AJ
     *
     * \param[in] AFrom point that must be appeared in \p AI and \p AJ corners
     * \param[in] AI    an end point
     * \param[in] AJ    another end point

     *
     * \return the id of common points
     */
    std::vector<int> findCommonPoints(const int AFrom,
                                      const int AI,
                                      const int AJ);
    /*------------------------------------------------------------------------*/
    /** \brief Starting from three corner points, it returns the free corner of
     *         the point that apperas in all of them and which is pointed to all
     *         of their origins
     *
     * \param[in]  ACorner1   first corner
     * \param[in]  ACorner2   second corner
     * \param[in]  ACorner3   third corner
     * \param[out] ACornerOut resulting corner
     *
     * \return true if we find the out corner, false otherwise
     */
    bool findCommmonLastCorner(const HexCorner& ACorner1,
                               const HexCorner& ACorner2,
                               const HexCorner& ACorner3,
                               HexCorner&       ACornerOut);
    /*------------------------------------------------------------------------*/
    /** \brief Looking at all the corners defined for origin point \p AOrigin,
     *         it returs the corner corresponding to (\p AOrigin, \p AI, \p AJ,
     *         \p AK) if it exist
     * \param[in]  AOrigin origin point
     * \param[in]  AI      an end point
     * \param[in]  AJ      another end point
     * \param[in]  AK      another end point
     * \param[out] AOut    the found hex corner if it exists
     *
     * \return true if we find the out corner, false otherwise
     */
    bool getCorner(const int AOrigin, const int AI, const int AJ,
                   const int AK, HexCorner& AOut);

    /*------------------------------------------------------------------------*/
    /** \brief Check if \p AC is a corner corresponding to (\p AOrigin, \p AI,
     *         \p AJ, \p AK)
     *
     * \param[in] AC      the found hex corner if it exists
     * \param[in] AOrigin origin point
     * \param[in] AI      an end point
     * \param[in] AJ      another end point
     * \param[in] AK      another end point
     *
     * \return true if they correspond, false otherwise
     */
    bool isCorner(const HexCorner& AC, const int AOrigin,
                  const int AI, const int AJ, const int AK);

    /*------------------------------------------------------------------------*/
    /** \brief Find among the corners in \p AIn those having \p AFrom as an
     *         adjacent point.
     *
     * \param[in]  AIn   the set of corners to be parsed
     * \param[in]  AFrom the point they must match
     * \param[out] AOut  corners of \p AIn having \p AFrom as end point
     *
     * \return the number of found compatible corners
     */
    int findFreeCorners(const std::vector<HexCorner>& AIn,
                        const int AFrom,
                        std::vector<HexCorner>& AOut);

    /*------------------------------------------------------------------------*/
    /** \brief Return the free corners of \p AOrigin with \p AI and \p AI as
     *         end points
     *
     * \param[in] AOrigin corner origin points
     * \param[in] AI      an end point
     * \param[in] AJ      another end point
     *
     * \return free corners of \p AOrigin with \p AI and \p AI as end points
     */
    std::vector<PointConnectionBuilder::HexCorner>
    findFreeCorners(const std::vector<int>& AOrigin,
                    const int AI, const int AJ);

    /*------------------------------------------------------------------------*/
    /** \brief Add a corner data to the point with index \p AIndex
     *
     * \param[in] AIndex  the point index we add a corner to
     * \param[in] AIndex1 the point index we go to in direction 1
     * \param[in] AV1     the vector in direction 1
     * \param[in] AIndex2 the point index we go to in direction 2
     * \param[in] AV2     the vector in direction 2
     * \param[in] AIndex3 the point index we go to in direction 3
     * \param[in] AV3     the vector in direction 3
     */
    void addCorner(const int AIndex,
                   const int AIndex1, const gmds::math::Vector3d& AV1,
                   const int AIndex2, const gmds::math::Vector3d& AV2,
                   const int AIndex3, const gmds::math::Vector3d& AV3);

    /*------------------------------------------------------------------------*/
    /** \brief Build all the corners associated to \p AOrigin by considering
     *         existing edges \p AEdges. This construction is based on geom.
     *         criteria
     *
     * \param[in] AOrigin origin point
     * \param[in] AEdges  existing edges connected to \p AOrigin
     */
    void buildCornersAsSolidAngles(const int AOrigin,
                                   std::vector<OrientedEdge>& AEdges);


    /*------------------------------------------------------------------------*/
    /** \brief Compute the mesh association for point AID. This computation is
     *         performed knowing how the point is classified (0,1,2 or 3) and
     *         the id of the curve
     *
     * \param[in] AID point id
     */
    void computeMeshAssociation(const int AID);

    /*------------------------------------------------------------------------*/
    /** \brief Return the region information (dim=3, id) of the region that
     *         contains \p APnt.
     * \param[in] APnt   the point we start from
     * \param[in] ATetID an initial tetrahedron we start from
     *
     * \return the info (3, id) that characterizes the tet containing APnt
     */
    gmds::Cell::Data getRegionContaining(const gmds::math::Point& APnt,
                                         const gmds::TCellID ATetID);
    /*------------------------------------------------------------------------*/
    /** \brief Return the face information (dim=2, id) of the boundary face that
     *         contains \p APnt.
     * \param[in] APnt    the point we start from
     * \param[in] ATetID  an initial tetrahedron we start from
     * \param[in] ASurfID the id of the surface \p APnt is supposed to be on
     *
     * \return the info (2, id) that characterizes the face containing APnt
     */
    gmds::Cell::Data getBoundaryFaceContaining(gmds::math::Point& APnt,
                                               const gmds::TCellID ATetID,
                                               const double ADistance,
                                               const int ASurfID);

    /*------------------------------------------------------------------------*/
    /** \brief Return the edge information (dim=2, id) of the boundary edge that
     *         contains \p APnt.
     * \param[in] APnt    the point we start from
     * \param[in] ATetID  an initial tetrahedron we start from
     * \param[in] ASurfID the id of the curve \p APnt is supposed to be on
     *
     * \return the info (1, id) that characterizes the edge containing APnt
     */
    gmds::Cell::Data getBoundaryEdgeContaining( gmds::math::Point& APnt,
                                               const gmds::TCellID ATetID,
                                               const double ADistance,
                                               const int ACurvID);
    /*------------------------------------------------------------------------*/
    /** \brief This function give the list of region ids such that generated
     *          points associated to those regions could be at a distance
     *          less or equal to \p AEps from  \p APnt
     *
     * \param[in] AP    the point we start from
     * \param[in] AT    the tetrahedron which should containt AP
     * \param[in] AEps  the distance we look for
     *
     * \return the ids of the tetrahedra that contains points at the specified
     *         distance
     */

    std::set<gmds::TCellID>
    getCloseRegionsFrom(const gmds::math::Point& AFromPnt,
                        const gmds::Region& AFromTet,
                        const double AEpsilon);
   static bool getOppositeRegion(const gmds::TCellID  ANodeID,
                           const gmds::Region&  AR,
                           gmds::Region&        AOut);
    static bool getOppositeFace(const gmds::TCellID  ANodeID,
                         const gmds::Region&  AR,
                         gmds::Face&        AOut);
    bool getOppositeBndFace(const gmds::TCellID ANodeID,
                            const gmds::Face&   AR,
                            gmds::Face&         AOut);
    void updateTetMesh();
    /*------------------------------------------------------------------------*/
    /** \brief Return the face of \p AR sharing nodes \p ANI and \p ANJ
     *         with face \p AFrom
     *
     * \param[in] AFrom  a face of \p AR whose \p ANI and \p ANJ are two nodes
     * \param[in] AR     a region
     * \param[in] ANI    a node id
     * \param[in] ANJ    a node id
     *
     * \return the expected face
     */

    static gmds::Face getFace(const gmds::Face& AFrom,
                       const gmds::Region& AR,
                       const gmds::TCellID ANI,
                       const gmds::TCellID ANJ);


    gmds::math::Vector3d getOutputNormal(gmds::Face& AFace, gmds::Region& ARegion);


    /*------------------------------------------------------------------------*/
    /** \brief Check if \p AI is in \p AV
     *
     * \param[in] AI an id
     * \param[in] AV a vector of ids
     *
     * \return true if \p AI is in \p AV, false otherwise
     */
    static bool isIn( gmds::TCellID AI, const std::vector<gmds::TCellID>& AV) ;

    static std::vector<gmds::TCellID> getFaces(const gmds::Node& ANI, const gmds::Node& ANJ) ;

    //The point is moved during this process
    gmds::Face closestFace(gmds::math::Point& AP, int ASurfID);
    gmds::Edge closestEdge(gmds::math::Point& AP,  int ACurvID);


    void writeInput();
    void writeHexes();
    void writeEdges(std::vector<std::vector<OrientedEdge> >& AInEdges,
                    const std::string& AFileName);
protected:
    /** Mesh we start from */
    gmds::Mesh* m_mesh;

    /** points used to build the mesh*/
    std::vector<gmds::math::Point> m_pnt;
    /** rotation associated to each point used to build the mesh*/
    const std::vector<gmds::math::Chart> m_chart;
    std::vector<gmds::Cell::Data> m_mesh_data;
    /** type associated to each point used to build the mesh*/
    const std::vector<int>               m_type;
    /** geometric classification of each point used to build the mesh*/
    const std::vector<int>               m_classification;
    /** Curve numbering of each point used to build the mesh*/
    const std::vector<int>               m_curve;
    /** surface numbering of each point used to build the mesh*/
    const std::vector<int>               m_surface;
    /** Defines a normal constraint alignment for boundary nodes
     * we work on */
     std::map<int, gmds::math::Vector3d> m_normal;

    /* vector given for each point if it is already used as a hex corner*/
    std::vector<int>               m_used;
    /** the generated hexahedral elements we wer mesh we work on */
    gmds::Mesh                   m_hexes;
    /** the expected spacing between points */
    double m_spacing;
    double m_dot_tolerance;
    /** we store computed edges in this variable */
    std::vector<std::vector<OrientedEdge> > m_edges;
    /** mapping from each point to the number of hex corner it must be in*/
    std::map<int, std::vector<HexCorner> > m_hc_mapping;
    std::map<int, gmds::Node> m_node_mapping;

    std::map<int,std::vector<int> > m_filter;

    bool m_with_debug_info;
    std::string m_output_dir;


};
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_FRAME3D_GRID_POINT_CONNECTOR_H_ */
/*----------------------------------------------------------------------------*/
