//
// Created by chenyt on 26/09/24.
//

#ifndef GMDS_NONCONFORMALHALFEDGE_H
#define GMDS_NONCONFORMALHALFEDGE_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include "gmds/medialaxis/MedialAxis2D.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API NonConformalHalfEdge
{
 private:
	// ID
	TCellID m_id;
	// Corresponding face
	Face m_face;
	// Corresponding conformal edges
	std::vector<Edge> m_conformal_edges;
	// First node
	Node m_first_node;
	// End node
	Node m_end_node;
	// Next half edge ID
	int m_next;
	// Opposite half edges IDs
	std::vector<int> m_opposite;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AID, AFace, AEdges
	 */
	NonConformalHalfEdge(TCellID AId, Face AFace, std::vector<Edge> AEdges);


	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~NonConformalHalfEdge()=default;

	/*-------------------------------------------------------------------------*/
	/** @brief Getters.
         *  @param
	 */
	TCellID  id(){return m_id;}
	Face face(){return m_face;}
	std::vector<Edge> edges(){return m_conformal_edges;}
	Node firstNode(){return m_first_node;}
	Node endNode(){return m_end_node;}
	int next(){return m_next;}
	std::vector<int> opposite(){return m_opposite;}

	/*-------------------------------------------------------------------------*/
	/** @brief Accessors.
         *  @param
	 */
	void next(int AN){m_next = AN;}
	void opposite(std::vector<int> AO){m_opposite = AO;}
	void firstNode(Node AN){m_first_node = AN;}
	void endNode(Node AN){m_end_node = AN;}

	/*-------------------------------------------------------------------------*/
	/** @brief Get the ordered nodes composing the half edge.
         *  @param
	 */
	std::vector<Node> getOrderedNodes();

};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_NONCONFORMALHALFEDGE_H
