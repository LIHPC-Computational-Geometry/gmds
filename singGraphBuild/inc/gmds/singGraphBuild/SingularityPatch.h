/*----------------------------------------------------------------------------*/
/*
 * SingularityPatch.h
 *
 *  Created on: June 4, 2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef SINGULARITYPATCH_H_
#define SINGULARITYPATCH_H_
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
class SingularityLine;
class SingularityPoint;
/*----------------------------------------------------------------------------*/
/* \class SingularityPatch
 * \brief Class that describes services and properties for singularity patches,
 *        that is a surface patch bounded by singularity lines connected in
 *        singularity points.
 */
class LIB_GMDS_SINGGRAPHBUILD_API SingularityPatch
{
 public:  

  /*------------------------------------------------------------------------*/
  /** \brief  Constructor
   */
  SingularityPatch();

  /*------------------------------------------------------------------------*/
  /** \brief  Destructor
   */
  virtual ~SingularityPatch();
  
  /*------------------------------------------------------------------------*/
  /** \brief Get the the list of points
   *
   * \param[OUT] APoint the list of singularity points
   * \return     true if the list is well ordered, false otherwise
   */
  bool getPoints(std::vector<SingularityPoint*>& APoints);
  
  /*------------------------------------------------------------------------*/
  /** \brief Get the the list of boundary lines
   *
   * \param[OUT] ALines the list of singularity lines
   * \return     true if the list is well ordered, false otherwise
   */
  bool getLines(std::vector<SingularityLine*>& ALines);

  /*------------------------------------------------------------------------*/
  /** \brief Add a singularity point to the boundary of (*this). If the point
   *         is already in the boundary loop it is not added.
   *
   * \param  APoint the point to connect with
   */
  void addPoint(SingularityPoint* APoint);

  /*------------------------------------------------------------------------*/
  /** \brief Add a singularity line to the boundary of (*this). If the line
   *         is already in the boundary loop it is not added.
   *
   * \param  ALine the line to connect with
   */
  void addLine(SingularityLine* ALine);

  
  /*------------------------------------------------------------------------*/
  /** \brief Check if the boundary of a patch is a single closed loop and 
   *         reorder the boundary lines and points if required.
   *
   *  \param return true if the boundary loop is valid, false otherwise
   */
  bool checkValidityAnReorder();

  /*------------------------------------------------------------------------*/
  /** \brief Get opposed line (line after two turn from AOpposedLine)
	*/
  SingularityLine *computeOpposedLine(SingularityLine *AOpposedLine);

 protected:

  /** flag indicating if the patch is safely oriented*/
  bool m_valid_flag;

  /** lines of the boundary */
  std::vector<SingularityLine*> m_lines;
  /** points of the boundary*/
  std::vector<SingularityPoint*> m_points;

  
};
/*-----------------------------------------------------------------*/
#endif /* SINGULARITYPATCH_H_ */
/*-----------------------------------------------------------------*/
