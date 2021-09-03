/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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
