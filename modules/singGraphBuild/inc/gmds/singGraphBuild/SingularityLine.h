/*----------------------------------------------------------------------------*/
/*
 * SingularityLine.h
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#ifndef SINGULARITYLINE_H_
#define SINGULARITYLINE_H_
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingularityPoint.h>
#include "LIB_GMDS_SING_GRAPH_BUILD_export.h"
/*----------------------------------------------------------------------------*/
class SingularityPatch;
class LIB_GMDS_SING_GRAPH_BUILD_API SingularityLine
{
 public:

  enum ESingularityGeomLineType {
    CURVE,	// sing. line located on a geometric curve or surf. It means 
    SURFACE, // frame field sing. pnt located on the surface
    VOLUME,		// frame field sing. pnt located in the volume
    GEOM_UNDEF	// frame field sing. pnt located in the volume
  };
  /*------------------------------------------------------------------------*/
  /** \brief  Default constructor.
   */
  SingularityLine(const  ESingularityGeomLineType AType = GEOM_UNDEF);

  /*------------------------------------------------------------------------*/
  /** \brief  Destructor.
   */
  virtual ~SingularityLine();
/*------------------------------------------------------------------------*/
  /** \brief  Indicates if the line is lyong on the mesh boundary
   */
  virtual bool isOnBoundary() const = 0;

  /*------------------------------------------------------------------------*/
  /** \brief  Smooth the line discretization (only meaningful in 3D volumes)
   */
  virtual void smooth(const int ANbStep=4) = 0;

  /*------------------------------------------------------------------------*/
  /** \brief  Returns the discrete length of the curve
   */
  double length() const;

  /*------------------------------------------------------------------------*/
  /** \brief  Returns the type of singularity line
   */
  ESingularityGeomLineType getType() const { return m_type; }

  /*------------------------------------------------------------------------*/
  /** \brief  Set the line number
   *
   * \param ANumber the new line number associated to (*this)
   */
  void setNumber(const int ANumber){m_number=ANumber;}

  /*------------------------------------------------------------------------*/
  /** \brief  Get access to the line number
   */
  int getNumber() const {return m_number;}
    
    /*------------------------------------------------------------------------*/
    /** \brief  Get the tangent to *this in AParam if it can compute it, i.e.
     *          *this is already discretized. Tangent is already given in
     *          direction of *this. It is also a unit vector!
     *
     * \param AParam the parameter we want to compute the tangent we want
     */
    gmds::math::Vector3d getTangent(const double AParam) const;
    /*------------------------------------------------------------------------*/
    /** \brief  Get the point location to *this in AParam if it can compute
     *          it, i.e. *this is already discretized. Tangent is already
     *          given indirection of *this.
     *
     * \param AParam the parameter we want to compute the tangent we want
     */
    gmds::math::Point getPoint(const double AParam) const;

  /*------------------------------------------------------------------------*/
  /** \brief Correct the line orientation. Considering the line discretization
   *         and its 2 end points, this operation invert the discretization
   *         to go from the 1st to the 2nd end point (if necessary).
   *
   * \return false if the operation is not performed (due to missing or wrong
   *         internal data). Ex: only one end point.
   */
  virtual bool healOrientation();
	/*------------------------------------------------------------------------*/
	/** \brief Add a end point, which is a singularity point. A line can be
		*         connected to 2 end points at most.
		*
		*  \param APnt a new end point of (*this)
		*/
	void addSlot(SingularityPoint::Slot *ASlot);
	void addSingularityPoint(SingularityPoint *ASing);

	/*------------------------------------------------------------------------*/
	/** \brief Remove the last added end point
		*
		*  \return the removed end point
		*/
	SingularityPoint::Slot *removeSlot();
  /*------------------------------------------------------------------------*/
  /** \brief Remove all connection to singularity points
   */
  void removeAllSingularityPoints();

  /*------------------------------------------------------------------------*/
  /** \brief Add a discretization point for the current line
   *
   *  \param APnt a geometric point
   */
  void addDiscretizationPoint(const gmds::math::Point& APnt);
  /*------------------------------------------------------------------------*/
  /** \brief Set the complete list  of discretization point
   *
   *  \param APnts new discretization list of points
   */
  void setDiscretizationPoints(const std::vector<gmds::math::Point >& APoints);
  /*------------------------------------------------------------------------*/
  /** \brief Return the discretization points of (*this)
   *
   *  \return a vector of geometric points
   */
  std::vector<gmds::math::Point >& getDiscretizationPoints();
  /*------------------------------------------------------------------------*/
  /** \brief Return the end points of (*this), i.e. 0, 1 or 2 singularity
   *         points
   *
   *  \return a vector of singularity points
   */
  std::vector<SingularityPoint* > getEndPoints() const;
  std::vector<SingularityPoint::Slot* >& getSlots();
  /*------------------------------------------------------------------------*/
  /** \brief Return the slot from the end point of (*this) at a specific location
	*         points
	*
	*  \param APointID index of slotLocation in the line discretization
	*  
	*  \return the slot at the given location
	*/
  SingularityPoint::Slot *getSlotAtLocation(const unsigned int &APointID);

	void setLinkedEdgeID(const int groupID) { m_linkedEdgeID = groupID;};
	int getLinkedEdgeID() { return m_linkedEdgeID;};

	void addPatch(SingularityPatch *patch) { m_patchs.push_back(patch);};
	std::vector<SingularityPatch *> getPatchs() { return m_patchs;};

 protected:
  ESingularityGeomLineType m_type;
  std::vector<SingularityPoint::Slot*> m_slots;
  std::vector<SingularityPatch *> m_patchs;
  std::vector<gmds::math::Point > m_discretization_points;
  int m_number;
  int m_linkedEdgeID = -1;

};
/*----------------------------------------------------------------------------*/
class VolumeSingularityLine : public SingularityLine
{
 public:

  /*------------------------------------------------------------------------*/
  /** \brief  Default constructor.
   */
 VolumeSingularityLine() :SingularityLine(SingularityLine::VOLUME){ ; }
 /*------------------------------------------------------------------------*/
  /** \brief Indicates if a line is on the boundary. Always false for volume
   *         lines.
   */
  bool isOnBoundary() const {return false;}

  virtual void smooth(const int ANbStep = 4);

};
/*----------------------------------------------------------------------------*/
class SurfaceSingularityLine : public SingularityLine
{
 public:

  /*------------------------------------------------------------------------*/
  /** \brief Default constructor.
   */
 SurfaceSingularityLine() :SingularityLine(SingularityLine::SURFACE){ ; }

  /*------------------------------------------------------------------------*/
  /** \brief Returns if the line is on the boundary. It is always true.
   */
  bool isOnBoundary() const { return true; }

  /*------------------------------------------------------------------------*/
  /** \brief Smooths the line to get a better shape.
   */
  virtual void smooth(const int ANbStep = 4);

  /*------------------------------------------------------------------------*/
  /** \brief Adds a face id into the list of traversed faces by (*this)
   *
   *  \param AID the face id we want to add
   */
  void addTraversedFace(const gmds::TCellID AID);
  
  /*------------------------------------------------------------------------*/
  /** \brief Returns a copy of the traversed faces id
   *
   *  \return the ids of the faces traversed by (*this)
   */
  std::vector<gmds::TCellID> getTraversedFaces() const;

  /*------------------------------------------------------------------------*/
  /** \brief Updates the full set of traversed faces id
   *
   *  \param AIDs the new set of face ids
   */
  void setTraversedFaces(const std::vector<gmds::TCellID>& AIDs);

  /*------------------------------------------------------------------------*/
  /** \brief Returns ture if AFaceID is the id of a traversed face. No specific
   *         ordering
   *
   * \param AFaceID the id of the face we want to check
   *
   *  \return true if the face is traversed, false otherwise
   */
  bool isTraversed(const gmds::TCellID AFaceID) const; 

  /*------------------------------------------------------------------------*/
  /** \brief Computes the intersection point between two surface singularity
   *         lines.
   *
   * \param[IN]  ALine the singularity line we want to compute the intersection 
   *             with.
   * \param[IN]  ARefPnt defines the area where we want to find the intersection
   * \param[IN]  ARefRadius defines the area where we want to find the 
   *             intersection
   * \param[OUT] APnt  the computed intersection point, if it exists.
   *
   *  \return true if there is an intersectio, false otherwise.
   */
  bool getIntersectionPoint( SurfaceSingularityLine* ALine, 
			     const gmds::math::Point& ARefPnt,
			     const double ARefRadius,
			     gmds::math::Point& APnt);
 protected:
  /** we store the ids of the faces that are traversed by the surface singularity
     line */
  std::vector<gmds::TCellID> m_traversed_faces_id;

};
/*----------------------------------------------------------------------------*/
class CurveSingularityLine : public SingularityLine
{
 public:
  /*------------------------------------------------------------------------*/
  /** \brief  Default constructor.
   */
 CurveSingularityLine() :SingularityLine(SingularityLine::CURVE){ ; }

  bool isOnBoundary() const { return true; }

  virtual void smooth(const int ANbStep = 4)
  { //wholly geometry-defined, thus impossible to smooth
    ; 
  }

  /*------------------------------------------------------------------------*/
  /** \brief Allows to set the edges that discretizes the curve
   */
  void setMeshEdges(const std::vector<gmds::Edge>& AEdges){
    m_ordered_mesh_edges = AEdges;
  }

  /*------------------------------------------------------------------------*/
  /** \brief Returns the set of edges that discretize the curve
   */
  std::vector <gmds::Edge>& getMeshEdges() {
    return m_ordered_mesh_edges;
  }

  /*------------------------------------------------------------------------*/
  /** \brief Correct the line orientation. Considering the line discretization
   *         and its 2 end points, this operation invert the discretization
   *         to go from the 1st to the 2nd end point (if necessary).
   *
   * \return false if the operation is not performed (due to missing or wrong
   *         internal data). Ex: only one end point.
   */
  virtual bool healOrientation();

 private:
  std::vector <gmds::Edge> m_ordered_mesh_edges;
};
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYLINE_H_ */
/*----------------------------------------------------------------------------*/
