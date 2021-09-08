/*----------------------------------------------------------------------------*/
#ifndef SH_FIELD_SOLVER_STRATEGY_ITF_H_
#define SH_FIELD_SOLVER_STRATEGY_ITF_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <map>
/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/utils/CommonTypes.h>
#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <gmds/math/SHarmonicL4.h>
#include <gmds/math/AxisAngleRotation.h>
#include "LIB_GMDS_FRAME_3D_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*---------------------------------------------------------------------------*/
// Frame File Headers
/*----------------------------------------------------------------------------*/
/** \class  FieldSoverStrategyItf
 *  \brief  Interface defining the operations that must be provided to the
 *          FieldGenerator cslass in order to build a 2D or 3D frame field
 *          (based on Spherical harmonics)
 */
class LIB_GMDS_FRAME_3D_API FieldSolverStrategyItf{

public:
    
    /*------------------------------------------------------------------------*/
    /** \brief Default constructor
     */
    FieldSolverStrategyItf():
    m_mesh(0),
    m_harmonic_field(0),
    m_chart_field(0),
    m_ordering(0),
    m_H0(0),
    m_H4(0),
    m_H8(0),
    m_surface_nodes(0),
    m_markNodeOnSurf(0),
    m_surface_mode(false),
    m_nb_unknowns(0),
    m_markNodeLocked(0),
    m_penalty(10)
    {;}
    
    /*------------------------------------------------------------------------*/
    /** \brief Default destructor
     */
    virtual ~FieldSolverStrategyItf(){};
    /*------------------------------------------------------------------------*/
    /** \brief Gives the mesh to work on
     *
     * \param[in] AMesh The mesh we work on
     */
    void setMesh(gmds::Mesh* AMesh){m_mesh = AMesh;}
    /*------------------------------------------------------------------------*/
    /** \brief Gives the Sph. Harmonic variable used to store SH in the mesh
     *
     * \param[in] AV the variable where SH are stored and updated
     */
    void setHarmonicVariable(gmds::Variable<gmds::math::SHarmonicL4>* AV){
        m_harmonic_field = AV;
    }
    /*------------------------------------------------------------------------*/
    /** \brief Gives the chart variable used to store orientations in the mesh
     *
     * \param[in] AV the variable where charts are stored and updated
     */
    void setChartVariable(gmds::Variable<gmds::math::Chart>* AV){
        m_chart_field = AV;
    }
    
    void setOrderingVariable(gmds::Variable<int>* AV){
        m_ordering = AV;
    }
    void setH0(std::map<gmds::TCellID, gmds::math::SHarmonicL4>* AMap){
        m_H0 = AMap;
    }
    void setH4(std::map<gmds::TCellID, gmds::math::SHarmonicL4>* AMap){
        m_H4 = AMap;
    }
    void setH8(std::map<gmds::TCellID, gmds::math::SHarmonicL4>* AMap){
        m_H8 = AMap;
    }
    void setNbUnknowns(const int ANb){m_nb_unknowns = ANb;}
    void setNbFreeNodes(const int ANb){m_nb_free_nodes = ANb;}
    void setMarkNodeLocked(const int AMark){m_markNodeLocked = AMark;}
    void setMarkNodeOnSurface(const int AMark){m_markNodeOnSurf = AMark;}
    void setPenaltyTerm(const double AD){m_penalty = AD;}
    void setSurfaceNodes(std::vector<gmds::Node>* AV){m_surface_nodes = AV;}
    void setSurfaceEdges(std::vector<gmds::Edge>* AE){m_surface_edges= AE;}
    void setSurfaceMode(const bool& ASurfMode){m_surface_mode = ASurfMode;}
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for solving the system
     */
    virtual void setX()=0;
    virtual void addSmoothingTerms()=0;
    virtual void addBoundaryConstraints()=0;
    virtual void addLockedTerms()=0;
    virtual void initializeAssembly()=0;
    virtual void finalizeAssembly()=0;
    virtual void addLocalOptimConstraints()=0;
    virtual void getFeasibleSolution()=0;
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for solving the system
     */
    virtual void solve()=0;

    virtual void clean()=0;

protected:
    Mesh* m_mesh;

    gmds::Variable<gmds::math::SHarmonicL4>* m_harmonic_field;
    gmds::Variable<gmds::math::Chart>* m_chart_field;
    gmds::Variable<int>* m_ordering;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4>* m_H0;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4>* m_H4;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4>* m_H8;
    std::vector<gmds::Node>* m_surface_nodes;
    std::vector<gmds::Edge>* m_surface_edges;
    bool m_surface_mode;
    int m_nb_unknowns;
    int m_nb_free_nodes;
    int m_markNodeLocked;
    int m_markNodeOnSurf;

    double m_penalty;
};
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_FIELD_SOLVER_STRATEGY_ITF_H_ */
/*----------------------------------------------------------------------------*/
