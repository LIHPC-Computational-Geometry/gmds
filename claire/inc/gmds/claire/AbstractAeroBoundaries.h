//
// Created by rochec on 23/03/2022.
//

#ifndef GMDS_ABSTRACTAEROBOUNDARIES_H
#define GMDS_ABSTRACTAEROBOUNDARIES_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*---------------------------------------------------------------------------*/
// Frame File Headers
/*----------------------------------------------------------------------------*/
/** \class  AbstractAeroBoundaries
 *  \brief  Classe abstraite pour les informations relatives aux frontières
 *  dans le cadre de l'aéro.
 */
class LIB_GMDS_CLAIRE_API AbstractAeroBoundaries{

 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;
	/*--------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit AbstractAeroBoundaries(Mesh *AMesh);
	/*--------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~AbstractAeroBoundaries();
	/*-------------------------------------------------------------------*/
	/** \brief Function to be called for initializing datas on the
	 * boundaries
	 */
	AbstractAeroBoundaries::STATUS execute();
	/*-------------------------------------------------------------------*/
	/** \brief Vérifie si un noeud est sur un bord.
	 */
	bool isBnd(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** \brief Vérifie si un noeud est sur la frontière amont ou non.
	 */
	bool isAmont(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** \brief Vérifie si un noeud est sur la paroi ou non.
	 */
	bool isParoi(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** \brief Récupère la marque sur les noeuds de bord.
	 */
	TInt getMarkBnd();
	/*-------------------------------------------------------------------*/
	/** \brief Récupère la marque sur les noeuds de la frontière amont.
	 */
	TInt getMarkAmont();
	/*-------------------------------------------------------------------*/
	/** \brief Récupère la marque sur les noeuds de la paroi.
	 */
	TInt getMarkParoi();
	/*-------------------------------------------------------------------*/
	/** \brief Vérifie si le/les objet/s sont immergés.
	 */
	bool isImmerged();
	/*-------------------------------------------------------------------*/
	/** \brief Retourne le nombre de bords.
	 */
	int getNbrBords();
	/*-------------------------------------------------------------------*/
	/** \brief Retourne la couleur de la frontière amont
	 */
	int getColorAmont();
	/*-------------------------------------------------------------------*/
	/** \brief Retourne la couleur du noeud.
	 */
	int getNodeColor(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Donne l'id d'un noeud "le plus à gauche" sur le bord de couleur
	 * color, correspondant en aéro au point d'arrêt.
	 */
	TCellID PointArret(int color);
	/*-------------------------------------------------------------------*/
	/** @brief Donne l'id du noeud le plus proche de n_id sur le bord
	 * de couleur color.
	 * @param color Couleur du bord regardé
	 * @param p Point auquel on cherche le noeud le plus proche
	 * @return Retourne l'id du noeud du bord color le plus proche du point
	 * p
	 */
	TCellID ClosestNodeOnBnd(int color, const math::Point& p);
	/*-------------------------------------------------------------------*/
	/** @brief Donne l'id d'un noeud sur le bord de couleur color.
	 * @param color Couleur du bord regardé
	 * @return Retourne l'id d'un noeud sur le bord color
	 */
	TCellID RandomNodeOnBnd(int color);
	/*-------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** \brief Marque les noeuds sur les bords dans la marque
	 * m_markBoundaryNodes
	 */
	virtual void MarkBoundariesNodes()=0;
	/*-------------------------------------------------------------------*/
	/** \brief Colorie les bords. Une couleur par bord connexe, 0 pour
	 * l'intérieur du domaine. Couleurs stockées dans m_var_color_bords.
	 */
	void ColoriageBordsConnexes();
	/*-------------------------------------------------------------------*/
	/** \brief Identifie de quelle couleur est le front Amont dans la
	 * variable m_var_color_bords.
	 */
	virtual void WhichColorIsAmont()=0;
	/*-------------------------------------------------------------------*/
	/** \brief Marque les noeuds sur la paroi dans la marque m_markNodesParoi
	 * et les noeuds sur la frontière amont dans la marque m_markNodesAmont.
	 */
	void MarkAmontAndParoiNodes();
	/*-------------------------------------------------------------------*/

 protected:
	/** pointer to mesh we work on */
	Mesh* m_mesh;
	/** Booléen qui vérifie si le maillage correspond à un/des objet/s immergé/s */
	bool m_isImmerged;
	/** Nombre de bords distincts */
	int m_nbrBords;
	/** Couleur qui correspond au bord extérieur */
	int m_color_Amont;
	/** Colorie chaque bord d'une couleur != */
	Variable<int>* m_var_color_bords ;
	/** Colorie chaque bord d'une couleur != */
	std::map<TCellID,int> m_map_color_bords;
	/** Vecteur des id des noeuds de bord */
	std::vector<TCellID> m_bnd_nodes_ids ;
	/** Marque les noeuds des bords */
	TInt m_markBoundaryNodes;
	/** Marque sur les noeuds de la paroi */
	TInt m_markNodesParoi;
	/** Marque sur les noeuds de la frontière amont */
	TInt m_markNodesAmont;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTAEROBOUNDARIES_H
