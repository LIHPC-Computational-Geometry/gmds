/*----------------------------------------------------------------------------*/
#ifndef SH_FIELD_GENERATOR_H_
#define SH_FIELD_GENERATOR_H_

/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <gmds/math/SHarmonicL4.h>
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include <gmds/frame3d/FieldSolverStrategyItf.h>
#include <gmds/frame3d/Params.h>
#include "LIB_GMDS_FRAME_3D_export.h"
/*---------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  FieldGenerator
 *  \brief  Computes a frame field onto a 3D mesh following a direct approach
 *          where a linear system is solved and frames are represented by
 *          Spherical Harmonics (SH)
 */
class LIB_GMDS_FRAME_3D_API FieldGenerator{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] ASolver the way the linear solver will be build and solve
     * \param[in] AMesh the mesh where we work on. Its a 3D mesh with regions
     * \param[in] AGParam parameters global to the full hex/quad meshing algo.
     * \param[in] AParam  parameters specific to the frame field generation
     * \param[in] ABM     parameters gathering boolean marks on mesh entities
     * \param[in] ASurfMode boolean indicating if we just work on the boundary
     */
    FieldGenerator(FieldSolverStrategyItf* ASolver,
                   gmds::Mesh* AMesh,
                   const ParamsGlobal& AGParam,
                   const ParamsFrameField& AParam,
                   const ParamsMark& ABM,
                   const bool ASurfMode=false);


    /*------------------------------------------------------------------------*/
    /** \brief Function to set a list of boundary nodes has having to be free
     *
     * @param ABndFreeNode variable which values true for a bnd node to be
     *                     free and false otherwise.
     */
    void relaxBoundary(Variable<int> *ABndFreeNode);

    void addHardConstraints(Variable<int> *AHasConstraint, Variable<math::Chart> *AConstraint);
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for generating the frame field
     */
    void execute();


    /*------------------------------------------------------------------------*/
    /** \brief Function returning the computed field energy as defined in
     *         RAY, SOKOLOF paper "On Smooth Frame Field"
     */
    double computeFieldEnergy();


    /*------------------------------------------------------------------------*/
    /** \brief Give access to the set of singular tets. A side effect of this
     *         operation is to update the region variable "sing_tet".
     *
     *  @return the list of singular tet ids (first) + for each tet its
     *          singularity type (second).
     */
    std::vector<std::pair<TCellID, int> > getSingularTetIds() const;

    /*------------------------------------------------------------------------*/
    /** \brief Returns the chart field build by this algorithm with 1 chart
     *         per mesh node.
     */
    Variable<gmds::math::Chart>* chartField(){
        return m_chart_field;
    }
protected:
    
    
    std::vector<gmds::Edge> getEdgesOnCurve(const gmds::Node& ANode) const;
    std::vector<gmds::Face> getFacesOnSurface(gmds::Edge& AEdge) const;
    gmds::Node  getNeighboorOn(const gmds::Node& AN, gmds::Edge& AE) const;
    
    
    
    virtual void buildFrameField() {}
    virtual void computeDistanceFields() {}
    virtual void initQuaternionsBySnapping() {}
    
    virtual void smoothAll(){}
    /*------------------------------------------------------------------------*/
    /** \brief Sort nodes to gather boundary nodes first then inner nodes
     */
    void sortNodes();
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute and store the normal to each boundary node. A single
     *         normal is assigned to nodes classified on surfaces, while a
     *         complete frame is computed for nodes classified on curves. Nodes
     *         classified on vertices are not taken into account.
     */
    void computeBoundaryData();

    void initHarmonicOnPoint  (const gmds::Node& ANode);
    void initHarmonicOnQuad  (const gmds::Node& ANode);
    void initHarmonicOnCurve  (const gmds::Node& ANode);
    void initHarmonicOnSurface(const gmds::Node& ANode);
    /*------------------------------------------------------------------------*/
    /** \brief Build the system to solve for iteration \p AI. Default value for
     *         \p AI is 0 meangin the intilialization stage. Previous solution
     *         is also provided if necessary.
     *
     * \param[in] APrevSolution a representation of the previous solution
     * \param[in] AI the iteration number for assembling the system
     */
    void buildSystem(const int AI=0);
    
  
    void writeSolution();
    gmds::math::Vector9d extractSolutionFromX(const gmds::Node& AN);

    
    /*------------------------------------------------------------------------*/
    /** \brief Initialize the boolean marks
     */
    void initMarks();
    
    /*------------------------------------------------------------------------*/
    /** \brief Clean the boolean marks
     */
    void cleanMarks();
    /*------------------------------------------------------------------------*/
    /** \brief Mark the boundary cells
     */
    void markBoundaryCells();

protected:
    
    /** the object used to solve the linear system */
    FieldSolverStrategyItf* m_solver;
    /** the mesh we work on */
    gmds::Mesh* m_mesh;
    
    ParamsGlobal m_global_params;
    ParamsFrameField m_params;
    ParamsMark m_bm;

    bool m_surface_mode;
    /*** Node mark for all the nodes classified on points and curves*/
    int m_markNodeFixed;
    /*** Node mark for all the nodes on the boundary*/
    int m_markNodeOnBnd;
    /*** Edge mark for all the edges used in the minimization pb*/
    int m_markEdgeSmooth;


    gmds::Variable<gmds::math::SHarmonicL4>* m_harmonic_field;
    gmds::Variable<gmds::math::Chart>*       m_chart_field;

    gmds::Variable<int>*       m_ordering;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4> m_H0;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4> m_H4;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4> m_H8;


    std::vector<gmds::Node> m_surface_nodes;
    std::vector<gmds::Edge> m_surface_edges;

    /*** number of inner, boundary nodes */
    int m_nb_inner_nodes; //nodes inside the volume
    int m_nb_fixed_nodes; //nodes on curves and pnts
    int m_nb_surface_nodes; //nodes on surface
    int m_nb_free_nodes; //inner+surface

    // All the edges, but those connecting two fixed nodes
    int m_nb_smooth_edges;
    // CURRENTLY UNUSED
    // int m_nb_strict_smooth_edges;

    /** after line correction, we get constraints to ensure*/
    bool m_with_hard_constraint;
    gmds::Variable<math::Chart> *m_hard_constraint;
    gmds::Variable<int> *m_has_hard_constraint;

    /** after line correction, we get constraints to ensure*/
    bool m_with_boundary_relax;
    gmds::Variable<int> *m_boundary_relax;

};
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_FIELD_GENERATOR_H_ */
/*----------------------------------------------------------------------------*/
