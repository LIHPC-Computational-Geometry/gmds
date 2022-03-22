//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE2D_H
#define GMDS_AEROPIPELINE2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/Params.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroPipeline2D
 *  \brief  Pipeline de génération de maillages 2D pour l'aéro.
 */
class LIB_GMDS_CLAIRE_API AeroPipeline2D: public AbstractAeroPipeline {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	AeroPipeline2D(ParamsAero Aparams);

	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual void execute();
	/*------------------------------------------------------------------------*/

 private:
	/*----------------------------------------------------------------------------*/
	/** @brief Lecture du fichier de maillage au format .vtk
	 */
	void LectureMaillage();
	/*----------------------------------------------------------------------------*/
	/** @brief Ecritures des fichiers de maillages triangulaires et quad au
	 * format .vtk
	 */
	void EcritureMaillage();
	/*----------------------------------------------------------------------------*/
	/** @brief Initialisation des marques sur les fronts
	 */
	void InitialisationFronts();
	/*----------------------------------------------------------------------------*/
	/** @brief Créé les sommets des blocs sur la paroi le maillage quad généré.
	 * Dans cette fonction, on récupère tous les noeuds du maillage triangulaire et
	 *
	 */
	void InitialisationMeshGen();
	/*----------------------------------------------------------------------------*/
	/** @brief Créé les sommets des blocs sur la paroi pour le maillage quad généré
	 */
	void DiscretisationBlocParoi();
	/*----------------------------------------------------------------------------*/
	/** @brief Génère une couche de noeuds du maillage
	 */
	void GenerationCouche(int couche_id, double dist);
	/*----------------------------------------------------------------------------*/
	/** @brief Donne le noeud de la couche i construit à partir du noeud n0 de la
	 * couche i-1. On suppose ici, pour l'instant, qu'il est unique.
	 */
	Node SuccessorNode(Node n0);
	/*----------------------------------------------------------------------------*/
	/** @brief Donne le noeud de la couche i construit à partir du noeud n0 de la
	 * couche i+1. On suppose ici, pour l'instant, qu'il est unique.
	 */
	Node AnteriorNode(Node n0);
	/*----------------------------------------------------------------------------*/
	/** @brief Donne l'id d'un noeud "le plus à gauche" sur le bord de couleur
	 * color, correspondant en aéro au point d'arrêt.
	 */
	TCellID PointArret(int color);
	/*----------------------------------------------------------------------------*/
	/** @brief Donne le vecteur des noeuds adjacents à n0 qui sont dans la couche i.
	 */
	std::vector<Node> AdjNodesInLayer(Node n0, int couche_i);
	/*----------------------------------------------------------------------------*/
	/** @brief Créé une face de type quad dans la couche avec les noeuds n0 de couche i,
	 * n1 son antécédant dans la couche i-1, n2 un noeud adjacent à n1 dans la couche i-1.
	 * Le dernier noeud n3 de la couche i est obtenu à l'aide de n2.
	 */
	void CreateQuadAndConnectivities(Node n0, Node n1, Node n2);
	/*----------------------------------------------------------------------------*/
	/** @brief Vérifie si une face de type quad est créée. Face correspondant
	 * aux noeuds n0 de couche i, n1 son antécédant dans la couche i-1, n2 un noeud
	 *  adjacent à n1 dans la couche i-1. Le dernier noeud n3 de la couche i est
	 *  obtenu à l'aide de n2.
	 */
	bool isQuadCreated(Node n0, Node n1, Node n2);
	/*----------------------------------------------------------------------------*/
	/** @brief Retourne la longueur d'un bord (périmètre) à partir d'un noeud au
	 * hasard
	 */
	double computeBoundaryLength(Node n0);
	/*----------------------------------------------------------------------------*/
	/** @brief Retourne la longueur d'un bord (périmètre) à partir de la couleur
	 * color du bord.
	 */
	double ComputeBoundaryLength(int color);
	/*----------------------------------------------------------------------------*/
	/** @brief Retourne un vecteur ordonné des noeuds du bord correspondant à la
	 * couleur color.
	 */
	std::vector<TCellID> BndNodesOrdered(int color);
	/*----------------------------------------------------------------------------*/
	/** @brief Créé les sommets des blocs sur le bord de couleur color pour le
	 * maillage quad généré.
	 */
	void DiscretisationBlocsBord(int color);
	/*----------------------------------------------------------------------------*/
 protected:
	/** mesh we work on */
	Mesh m_m;
	/** mesh quad generated */
	Mesh m_mGen;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE2D_H
