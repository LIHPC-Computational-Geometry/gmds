/*----------------------------------------------------------------------------*/
#ifndef PGP_COMPUTING_H_
#define PGP_COMPUTING_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
#include <set>
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <gmds/math/SHarmonicL4.h>
/*----------------------------------------------------------------------------*/
#include <gmds/frame3d/Params.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** @class  PGPComputing
 *  @brief  Starting from a frame field, a set of points and possible connection
 *          between them are extracted. The algorithms relies on working in
 *          a parametrization that fits the input frame field as best as possible.
 */
class EXPORT_GMDS PGPComputing{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] AMesh    the mesh where we work on
     * \param[in] AMarks   Boolean marks
     * \param[in] ASpacing the expected spacing between generated points.
     * \param[in] ACurl    used curl value for the algorithm
     */
    PGPComputing(gmds::Mesh* AMesh,
                 const ParamsGlobal& AParamGl,
                 const ParamsMark& AMarks,
                 const double ASpacing=1.0,
                 const double ACurl=0.35);


    std::map<gmds::TCellID, gmds::math::Vector3d >& getUI() {return m_Ui;}
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for generating the frame field
     */
    void execute();


protected:

    void computeCurlCorrection();
    /*------------------------------------------------------------------------*/
    /** \brief build the local id map to have a continous numbering of the nodes
     *        involved in this algorithm
     */
    void buildLocalIDs();

    /*------------------------------------------------------------------------*/
    /** \brief Nodes located on the boundary are constrained by normal
     *         preserving
     */
    void setBoundaryConstraint();

    /*------------------------------------------------------------------------*/
    /** \brief Init the least-square system to be solved
     */
    void initSystem();
    /*------------------------------------------------------------------------*/
    /** \brief Build the least-square system to be solved
     */
    void buildSystem();
    /*------------------------------------------------------------------------*/
    /** \brief Solve the least-square system to be solved
     */
    void solveSystem();
    /*------------------------------------------------------------------------*/
    /** \brief Extract the parametrization solution after the system resolution
     */
    void getUiSolution();
    /*------------------------------------------------------------------------*/
    /** \brief Clean up the memory used by the least-square system solver
     */
    void cleanSystem();

    /*------------------------------------------------------------------------*/
    /** \brief Compute the the g_ij terme for the oriented edge defined by \p
     *         AEdge and going from node \p AFrom to node \p ATo
     *
     * \param[in] AEdge an  edge
     * \param[in] AFrom the node we go from
     * \param[in] ATo   the node we go to
     *
     * \param the value Gij computed for this oriented edge
     */
    gmds::math::Vector3d computeGijWithCurl(gmds::Edge& AEdge,
                                            gmds::Node& AFrom,
                                            gmds::Node& ATo);

    /*------------------------------------------------------------------------*/
    /** \brief Structure used to orient edges for Rij computation
     */
    struct OrientedEdge{
        /** mesh edge this is relative too*/
        gmds::Edge edge;
        /** Starting point of the oriented edge*/
        gmds::Node first;
        /** End point of the oriented edge*/
        gmds::Node second;
        /*--------------------------------------------------------------*/
        /** \brief Default constructor (for STL container purpose only)
         */
        OrientedEdge();
        /*--------------------------------------------------------------*/
        /** \brief Constructor from an edge and two nodes
         */
        OrientedEdge(gmds::Edge& e,
                     gmds::Node& f,
                     gmds::Node& s);

        /*--------------------------------------------------------------*/
        /** \brief Constructor from an edge only
         */
        OrientedEdge(gmds::Edge& e);

        /*--------------------------------------------------------------*/
        /** \brief Returns true if this->first and this->second are resp.
         *         the first and second node of this->edge, false otherwise
         */
        bool isWellOriented();
    };

    /*------------------------------------------------------------------------*/
    /** \brief Build the oriented edges of face \p AF in such a way that the
     *         edge i connects node i to i+1 (modulo the number of corners)
     *
     * \param[in] AF a face
     *
     * \return an ordered set of oriented edges
     */
    std::vector<OrientedEdge> getOrientedEdge(const gmds::Face& AF);


    /*------------------------------------------------------------------------*/
    /** \brief Compute the g_ij value for oriented edge \p AEdge
     *
     * \param[in] AEdge an oriented edge
     *
     * \return a 3D vector
     */
    gmds::math::Vector3d computeGij(OrientedEdge& AEdge) ;

    /*------------------------------------------------------------------------*/
    /** \brief Compute the c_ij corrective term for oriented edge \p AEdge
     *
     * \param[in] AEdge an oriented edge
     *
     * \return a 3D vector
     */
    gmds::math::Vector3d computeCij(OrientedEdge& AEdge) ;

    /*------------------------------------------------------------------------*/
    /** \brief Compute the Rij mapping transporting the frame defined in \p AI
               into the frame defined in \p AJ
     *
     * \param[in]  AI a first node id
     * \param[in]  AJ a second node id
     * \return a mapping transporting the frame defined in \p AI into the
     *         frame defined in \p AJ
     */
    gmds::math::Chart::Mapping getRij(const gmds::TCellID AI,
                                      const gmds::TCellID AJ) const;
    /*------------------------------------------------------------------------*/
    /**  \brief Compute the Rij mapping transporting the frame defined in \p AI
     *          into the frame defined in \p AJ
     *
     * \param[in]  AI a first node
     * \param[in]  AJ a second node
     * \return a mapping transporting the frame defined in \p AI into the
     *         frame defined in \p AJ
     */
    gmds::math::Chart::Mapping getRij(const gmds::Node& AI,
                                      const gmds::Node& AJ) const;

    /*------------------------------------------------------------------------*/
    /**  \brief This method align the Tij value \p AToBeAligned onto the Tij
     *          value \p ARef, knowing that a geometrical deviation must be
     *          taken into account
     *
     * \param[in/out] AToBeAligned   a Tij value to align
     * \param[in]     ARef           the reference Tij value
     * \param[in]     AGeomDeviation the known geometric deviation between the
     *                               two first parameters
     */
    void alignTij(gmds::math::Vector3d& AToBeAligned,
                  const gmds::math::Vector3d& ARef,
                  const gmds::math::Vector3d& AGeomDeviation);

protected:
    
    /** the mesh we work on */
    gmds::Mesh* m_mesh;
    /** global parameter including parameters for this algorithm*/
    ParamsGlobal m_param_gl;

    /** rotation field previously computed for each node. */
    gmds::Variable<gmds::math::AxisAngleRotation>* m_rotation_field;


    /** Boolean marks */
    ParamsMark m_bm;

    /** the relative spacing distance between generated points */
    double m_spacing;

    double m_curl;

    /** nb unknowns of the least-square system to be solved */
    int m_nb_unknowns;


    /** Give a local id for the nodes we work on in this algorithm */
    std::map<gmds::TCellID, int> m_id;

    /** Give a local id for the edges we work on in this algorithm */
    std::map<gmds::TCellID, int> m_edge_id;

    std::map<gmds::TCellID, gmds::math::Vector3d > m_Ui;
    /** correction factor along each edge. The global id of the edge is
     * used as the key */
    std::map<gmds::TCellID, gmds::math::Vector3d > m_corr;

    /** For a boundary node indicate which of its chart vector is
     * aligned with the surface normal N. A vector (0,1,0) indicates
     * that the second chart vector is aligned with N*/
    std::map<gmds::TCellID, gmds::math::Vector3d > m_bnd_constraint;

};

    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_EDGE_GENERATOR_H_ */
/*----------------------------------------------------------------------------*/
