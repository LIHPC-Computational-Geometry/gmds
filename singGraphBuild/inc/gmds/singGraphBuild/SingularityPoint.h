/*----------------------------------------------------------------------------*/
/*
 * SingularityPoint.h
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#ifndef SINGULARITYPOINT_H_
#define SINGULARITYPOINT_H_
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
class SingularityLine;
/*----------------------------------------------------------------------------*/
/* \class SingularityPoint
 * \brief Class that describes services and properties shared by
 *        all types of singularity points.
 */
class EXPORT_GMDS SingularityPoint
{
 public:  
  /*-------------------------------------------------------------------------*/
  /* \brief A Slot define local connection data between a singularity
   *        point and singularity lines
   *
   * More specifically, for each slot, we store
   * - if a singularity line is already connected to this slot
   * - the starting point of the line for this singularity point
   * - the direction of the line for this singularity point
   * - the cell the line start from *this. It can be an edge or a node
   * - the line itself (if it has been launched previously)
   * - the singularity point itself
   */
  struct Slot{
    bool isLaunched;
    gmds::math::Point location;
    gmds::math::Vector3d direction;
    gmds::TCellID starting_cell_id; //not relevant for geom sing. point
    int  starting_cell_dim;         //not relevant for geom sing. point
    SingularityLine* line;
    gmds::math::Vector3d line_direction;
    SingularityPoint* from_point;
    bool isOnSurface;
    bool isFreeze; //means we can not disconnect a curve for it
    double lineDeviation; // the deviation along the singularity line connected to the slot
    Slot();
  };
  /*------------------------------------------------------------------------*/
  /* \brief Type of a singularity point
   */
  enum ESingularityType {
    FIELD, // singularity of the frame field 
    GEOM   // singularity in the geometry (thus regular frame field around it)
  };
  /*------------------------------------------------------------------------*/
  /* \brief Geometric type of a singularity point. 
   */
  enum ESingularityGeomType {
    VERTEX,   // Sing point located onto a geometric vertex
    CURVE,    // Sing point located on a geometric curve or surf. 
    SURFACE,  // Frame field sing. pnt located on the surface
    VOLUME,   // Frame field sing. pnt located in the volume
    GEOM_UNDEF// ????
  };

  /*------------------------------------------------------------------------*/
  /** \brief  Constructor
   *
   * \param AOwner the background mesh the singularity structure lives on
   * \param AType  the singularity type
   * \param AGType the geometric singularity type
   */
  SingularityPoint(gmds::Mesh* AOwner,
		   const ESingularityType AType = FIELD,
		   const ESingularityGeomType AGType = GEOM_UNDEF);
  /*------------------------------------------------------------------------*/
  /** \brief  Constructor
   *
   * \param AOwner the background mesh the singularity structure lives on
   * \param AX     the X coordinate of the singularity location
   * \param AY     the Y coordinate of the singularity location
   * \param AZ     the Z coordinate of the singularity location
   * \param AType  the singularity type
   * \param AGType the geometric singularity type
   */
  SingularityPoint(gmds::Mesh* AOwner,
		   const double AX, const double AY, const double AZ,
		   const ESingularityType AType = FIELD,
		   const ESingularityGeomType AGType = GEOM_UNDEF);

  /*------------------------------------------------------------------------*/
  /** \brief  Destructor
   */
  virtual ~SingularityPoint();

  /*------------------------------------------------------------------------*/
  /** \brief  Return the type of the singularity point
   */
  virtual ESingularityType getType() const {
    return m_type;
  } 

  /*------------------------------------------------------------------------*/
  /** \brief  Return the  geometric type of the singularity point
   */
  virtual ESingularityGeomType getGeomType() const 
  {
    return m_geom_type;
  }
  /*------------------------------------------------------------------------*/
  /** \brief  Change the location of the singularity point
   *
   * \param AX     the X coordinate of the singularity location
   * \param AY     the Y coordinate of the singularity location
   * \param AZ     the Z coordinate of the singularity location
   */
  void setXYZ(const double AX, const double AY, const double AZ);

  /*------------------------------------------------------------------------*/
  /** \brief  Change the location of the singularity point
   *
   * \param APnt the new location
   */
  void setLocation(const gmds::math::Point& APnt);
  
  /*------------------------------------------------------------------------*/
  /** \brief Get the location of the singularity point
   *
   * \return a location point
   */
  gmds::math::Point getLocation() const;

  /*------------------------------------------------------------------------*/
  /** \brief Connect a singularity line to (*this). In this case, it means 
   *         finding an open slot along direction AVec.
   *
   * \param  ALine the line to connect with
   * \param  AVec  the direction the line arrives on (*this)
   * \return true if we are able to connect, false otherwise (due to a bad
   *         direction (AVec) or no available slot
   */
  virtual bool connectLine(SingularityLine* ALine,
		       const gmds::math::Vector3d& AVec);
  
  /*------------------------------------------------------------------------*/
  /** \brief Connect a singularity line to (*this). Warning, when we add a line,
   *         the connection is well performed only if the line discretization is
   *         already known.
   *
   * \param  ALine the line to connect with
   */
  virtual void addLine(SingularityLine* ALine);

  /*------------------------------------------------------------------------*/
  /** \brief Creates a new slot
   */
  Slot* newSlot(const gmds::math::Point& APnt,
		const gmds::math::Vector3d& AVec,
		const gmds::TCellID& ACellID, 
		const int ACellDim,
		const bool AIsOnSurf,
		SingularityLine* ALine = 0,
		const gmds::math::Vector3d ALineDirection = gmds::math::Vector3d(1,0,0));

  /*------------------------------------------------------------------------*/
  /** \brief Creates a new geometric slot.  A geom slot is associated to a curve
   *         that is geometrically defined
   */
  Slot* newGeomSlot(const gmds::math::Point& APnt, 
		    SingularityLine* ALine = 0);

  /*------------------------------------------------------------------------*/
  /** \brief Removes all the slot of (*this) and so disconnect (*this) from all
   *         the connected curves. The curve->point connection is not updated
   */
  void clearSlots();

  /*------------------------------------------------------------------------*/
  /** \brief Get the list of all slots
   */
  std::vector<Slot*>& getSlots() { return m_slots; }
  /*------------------------------------------------------------------------*/
  /** \brief Get the list of all the lines connected to (*this)
   */
  std::vector<SingularityLine*> getLines();

  /*------------------------------------------------------------------------*/
  /** \brief Returns the singularity index
   */
  int index()const { return m_slots.size(); }
  /*------------------------------------------------------------------------*/
  /** \brief Returns the number of mesh cells the singularity is associated to.
   */
  int getNbMeshCells() const  {return m_cell_id_owner.size(); }

 /*------------------------------------------------------------------------*/
  /** \brief Return the mesh cell of type T the singularity point is 
   *         associated to. It can be 0, 1 or more depending on the type of 
   *	     singularity.
   *
   * \return a collection of T-type  cells
   */
  template<class T> std::vector<T> getMesh();

 /*------------------------------------------------------------------------*/
  /** \brief Returns the next line on the surface considering the singularity
   *         is not inside a volume and that ALine is on a curve or surface
   *
   * \param  ALine a line classified onto a curve or a surface
   * \return the next singularity line incident to *this and on classified on
   *         a curve or surface 
   */
  SingularityLine* nextLine(SingularityLine* ALine);

  void addMeshNode(const gmds::Node& ANode);
  void addMeshEdge(const gmds::Edge& AEdge);
  void addMeshFace(gmds::Face& AFace);
  void addMeshRegion(gmds::Region& ARegion);
 protected:

  /** background mesh the point is assigned to */
  gmds::Mesh* m_mesh; 
  
  /** Singularity type */
  ESingularityType m_type;
  
  /** Singularity geometric type */
  ESingularityGeomType m_geom_type;
      
  /** degree of the singularity, i.e. its number of slots */
  int m_degree;
  
  /** Collection of slots (free or connected)*/
  std::vector<Slot*> m_slots;
  
  /** geometric location */
  gmds::math::Point m_location;

  /** mesh cells that contains the singularity (it can be any type of mesh cells
   *   (node, edge, face or region) */
  /** mesh cell ids */
  std::vector<gmds::TCellID> m_cell_id_owner;
  /** mesh cell dimension */
  std::vector<int> m_dim_owner;
};
/*-----------------------------------------------------------------------------*/
/* \class SurfaceSingularityPoint
 *
 * \brief Specialization of singularity points for points living in faces. It
 *        can be a singularity point of the frame field or the intersection of
 *        two singularity lines in a stable area.
 */
class EXPORT_GMDS SurfaceSingularityPoint : public SingularityPoint {
 public:

  /*------------------------------------------------------------------------*/
  /** \brief Constructor
   */
 SurfaceSingularityPoint(
			 gmds::Mesh* AOwner,
			 const double AX, const double AY, const double AZ,
			 const ESingularityType AType=FIELD) 
   :SingularityPoint(AOwner,AX,AY,AZ,AType,SURFACE){; }
  /*------------------------------------------------------------------------*/
  /** \brief Constructor
   */
 SurfaceSingularityPoint(gmds::Mesh* AOwner, const ESingularityType AType = FIELD)
   :SingularityPoint(AOwner,AType,SURFACE)
    {;}

  /*------------------------------------------------------------------------*/
  /** \brief Returns the mesh face that contains this singularity
   *
   * \return a mesh face
   */
  gmds::Face getMeshFace();
 /*------------------------------------------------------------------------*/
  /** \brief Add a line incoming from the adjacent volume
   *
   * \param ALine the line to be addeds
   */
  bool addLineFromVolume(SingularityLine* ALine);

};
/*-----------------------------------------------------------------*/
/* \class CurveSingularityPoint
 * \brief This singularity point lies onto a geometrical curve
 *	  As a consequence, its slots corresponds to the geometrical
 *        curve it lies on in 2 directions + 1 or 2 line inside the
 *        adjacent surfaces.
 */
class EXPORT_GMDS CurveSingularityPoint : public SingularityPoint {
 public:

 CurveSingularityPoint(gmds::Mesh* AOwner, 
		       const double AX, const double AY, const double AZ,
		       const ESingularityType AType = GEOM)
   :SingularityPoint(AOwner, AX, AY, AZ, AType, CURVE){;}

 CurveSingularityPoint(gmds::Mesh* AOwner, const ESingularityType AType = GEOM)
   :SingularityPoint(AOwner,AType, CURVE){;}


};
/*-----------------------------------------------------------------*/
/* \class VertexSingularityPoint
 * \brief This singularity point corresponds to a geometrical
 *		  vertex V.
 *		  As a consequence, its slots corresponds to the geometrical
 *        curves incident to V or to inner volume sing. line
 */
class EXPORT_GMDS VertexSingularityPoint : public SingularityPoint {
 public:

 VertexSingularityPoint(gmds::Mesh* AOwner, 
			const double AX, const double AY, const double AZ,
			const ESingularityType AType = GEOM)
   :SingularityPoint(AOwner,AX, AY, AZ, AType, VERTEX){;}

 VertexSingularityPoint(gmds::Mesh* AOwner, 
			const ESingularityType AType = GEOM)
   :SingularityPoint(AOwner, AType, VERTEX){;}


  gmds::Node getMeshNode();
};
/*-----------------------------------------------------------------*/
class EXPORT_GMDS VolumeSingularityPoint : public SingularityPoint {
 public:
 VolumeSingularityPoint(gmds::Mesh* AOwner) 
   :SingularityPoint(AOwner,FIELD, VOLUME){ ; }

  virtual bool addLine(SingularityLine* ALine,
		       const gmds::Face& AIncomingFace);

  gmds::Region getMeshRegion();
  std::vector<gmds::Region> getMeshRegions();
 private:

  //list of tetra in the cluster of sing around *this
  std::vector<int> m_clusterTetID;

};
/*-----------------------------------------------------------------*/
template<> EXPORT_GMDS std::vector<gmds::Node>   SingularityPoint::getMesh<gmds::Node>();
template<> EXPORT_GMDS std::vector<gmds::Edge>   SingularityPoint::getMesh<gmds::Edge>();
template<> EXPORT_GMDS std::vector<gmds::Face>   SingularityPoint::getMesh<gmds::Face>();
template<> EXPORT_GMDS std::vector<gmds::Region> SingularityPoint::getMesh<gmds::Region>();
/*-----------------------------------------------------------------*/
#endif /* SINGULARITYPOINT_H_ */
/*-----------------------------------------------------------------*/
