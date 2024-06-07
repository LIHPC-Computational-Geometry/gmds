/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    InterfaceNodesPosSmoothVF.cpp
 *  \author  legoff
 *  \date    09/30/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
#include <random>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <glpk.h>

#include <GCoptimization.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/VectorDyn.h>

#include <KM/Utils/Graph.h>
#include <KM/DS/Mesh.h>
#include <KM/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
#include <ELG3D/DATACMPNT/Parameters.h>
#include <ELG3D/ALGOCMPNT/Tools.h>

#include <ortools/linear_solver/linear_solver.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace elg3d;
namespace milp {

    void
    InterfaceNodesPosSmoothVF_assignORTOOLS_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                           const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                           const FracPres *Afp_source,
                                           MaterialAssignment *Ama_target,
                                           const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                           const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                           const kmds::Variable<bool> *AVarMixedCells_source,
                                           const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                           const kmds::Graph *AGraph,
                                           const int ANbSubPixels)
    {
	    // =======================================
	    // Initialisation des variables, avec nbVert= nb de sommets dans un graph , nbMat = nb de matériaux, nbCells_target= nb de cellules cibles
	    const int nbVert = AGraph->getNbVec();              // récupère nb de sommets dans le graph Agraph
	    const int nbMat = Afp_source->getNbMaterials();     // récupère nb de matériel l'objet Afp_source

	    const int nbCells_target = ACellsIDs_target->getNbElems();     // on récupère le nombre de cellules cibles dans l'objet ACellsIDs_target

	    // ici on va créer une map pour les pixels purs.Le but est d'associer chaque pixel à un matériau pur. Grâce à cette étape on peut savoir quels pixels sont purs et à quels matériaux ils sont associés

	    std::map<kmds::TCellID, int> vert2puremat;     // on va créer une map qui associe chaque sommet(kmds::TCellID) à un entier représentant le matériau pur.
	    for (int i = 0; i < nbCells_target; i++) {           // boucle qui parcours chaque cellule cible
		    kmds::TCellID cid = ACellsIDs_target->get(i);     // récupère l'id de la cellule cible actuelle

		    kmds::TCellID vert = (*AVarCells2vertices)[cid];                 // récupère le sommet associé à cette cellule
		    if (vert != kmds::NullID) {                                      // vérifie si le sommet est valide
			    kmds::TCellID cid_old = (*AVarSubCells2oldCells)[cid];        // récupère l'ancienne cellule associée
			    if (!Afp_source->isMixedCell(cid_old)) {                      // vérifie si cette ancienne cellule n'est pas mixte
				    for (int imat = 0; imat < nbMat; imat++) {                 // parcourt chaque matériau
					    if (Afp_source->getFracPres(imat, cid_old) == 1.) {     // vérifie si le matériau imat a une fraction de cellule cid_old. cela signifie que imat est le matériau pur de cette cellule
						    vert2puremat.emplace(
						       vert,
						       imat);     // si la condition prcdt est remplie, ajoute l'association entre le sommet vert et le matériau imat à la map vert2puremat
						    break;
					    }
				    }
			    }
		    }
	    }

	    // calcul du nombre de cellules mixtes dans le graphe. cellules mixtes = cellules qui possèdeent plrs matériaux
	    int nbMixedCells_source = 0;     // initialise  un compteur pour le nombre de cellules mixtes
	    Kokkos::parallel_reduce(
	       ACellsIDs_source->getNbElems(),     // utilise une réduction en parallèle pour compter les cellules / récupère le nombre de cellules sources
	       KOKKOS_LAMBDA(const int i, int &sum) {               // la fonction lambda qui prend en indice i et une réfécrence à sum
		       kmds::TCellID cid = ACellsIDs_source->get(i);     // récupère l'id de la cellule source actuelle
		       if ((*AVarMixedCells_source)[cid]) {              // incrémente sum si la cellule est mixte
			       sum++;
		       }
	       },
	       nbMixedCells_source);
/*

	    // code convertit
	    //  =======================================
	    //  ORTOOLS
	    //  =======================================
	    // Initialisation de GLPK/on cherche aussi à minimiser la fonction objective
	    operations_research::MPSolver solver("InterfaceNodesPosSmoothVF", operations_research::MPSolver::SCIP_MIXED_INTEGER_PROGRAMMING);
	   /* if (!solver) {
		    std::cerr << "echec dans la création du solver." << std::endl;
		    return;
	    }
*/
	    // compute problem size
	    // Détermination du nb total de variables et de contraintes
	    int nbRows = 0;     // initialise le compteur du nb de lignes à 0
	    int nbCols = 0;     // initialise le compteur de nb de colonnes à 0
	    int nnz = 0;        // initialise le compteur du nb d'éléments non nuls dans la matrice des contraintes à 0

	    // pixel color variables/variables de couleur des pixels
	    nbCols += nbVert * nbMat;     // ajoute au compteur de colonnes le nb tot de variables de couleur des pixels, soit nbVert*nbMat

	    // color unicity per pixel/unicité de la couleur par pixel: 1 pixel ne peut avoir qu'une seule couleur attribuée
	    nbRows += nbVert;          // ajout au compteur  de ligne une contrainte pour chaque pixel, -> chaque pixel n'a qu'une seule couleure
	    nnz += nbVert * nbMat;     // ajoute au compteur d'éléments != 0 le produit du nb de pixels par le nb de matériaux

	    // vf per coarse mixed cell	/ calcul du pourcentage de matériaux  par cellule mixte: la somme des fractions volumiques des matériaux dans une seule cellule est égale à un nombre spécifique.
	    nbRows += nbMixedCells_source
	              * nbMat;     // ajoute au compteur lignes des contraintes pour les fractions volumiques de chaque cellule mixte et chaque matériaux
	    nnz += nbMixedCells_source * nbMat * ANbSubPixels;     // ajoute au compteur d'éléments !=0 les contributions de chaque cellule mixte et de chaque matériau multiplié par le nombre de sius-pixels

	    // neighboring similarity/ cette contrainte assure que les couleurs des pixels voisins sont les mêmes
	    nbRows += nbVert * nbMat;         // ajoute au compteur de ligne des contraintes pour la similarité de voisinge entre les pixels.
	    nbCols += 2 * nbVert * nbMat;     // ajoute au compteur des colonnes des variables supplémentaires pour les contraintes de voisinage
	    nnz += 3 * nbVert * nbMat
	           + nbMat * AGraph->getNbEdges();     // ajoute au compteur d'éléments des contraintes de similarité de voisinage et des arêtes du graphe des pixels

	    // ajouts des contraintes des variables dans GLPK
	    // pure pixels
	    nbRows += vert2puremat.size();     // ajoute au compteur de lignes un contrinte pour chaque pixel pur( où la couleur est fixée à un matériau spécifique)
	    nnz += vert2puremat.size();     // ajoute au compteur d'éléments non nuls les contributions de chaque pixel pur.


/*
	    //===================================================================================================================================================================================================================================================
	    // enforce constraints
	    glp_add_rows(lp, nbRows);     //
	    glp_add_cols(lp, nbCols);

	    ia = (int *) calloc(1 + nnz, sizeof(int));
	    ja = (int *) calloc(1 + nnz, sizeof(int));
	    ar = (double *) calloc(1 + nnz, sizeof(double));

	    int irow = 1;     // start at 1, not 0 !
	    int icol = 1;     // start at 1, not 0 !
	    int innz = 1;     // start at 1, not 0 !

	    // pixel color variables/contraintes de couleur
	    // represented as booleans
	    for (int p = 0; p < nbVert; p++) {     // boule parcourue pour chaque pixel et matériau
		    for (int imat = 0; imat < nbMat; imat++) {
			    std::string col_name =
			       "pixelcolor_" + std::to_string(p) + "_" + std::to_string(imat);     // crée un nom unique pour chaque variable de couleur de pixel
			    glp_set_col_name(lp, icol, col_name.c_str());                          // définit le nom de la variable dans le problème GLPK
			    glp_set_col_kind(lp, icol, GLP_BV);                                    // définit la variable comme étant binaire
			    // glp_set_col_bnds(lp, icol, GLP_DB, 0, 1); // not necessary because set by GLP_BV

			    glp_set_obj_coef(
			       lp, icol,
			       0.);     // définit le coefficient de l'objectif pour cette variable à 0 (cad qu'elle n'affecte pas directement l'objectif à minimiser)

			    icol++;     // incrémente le compteur de colonnes

			    // mettre le code convertit ici
			    // declarer un tableau ou je mets les valeurs et les variables sans déranger le code
		    }
	    }
	    */
	    // code convertit
	    //  =======================================
	    //  ORTOOLS
	    //  =======================================
	    std::vector<operations_research::MPVariable*> variables;
	    for (int p = 0; p < nbVert; p++) {
		    for (int imat = 0; imat < nbMat; imat++) {
			    std::string col_name = "pixelColor_" + std::to_string(p) + "_" + std::to_string(imat);
			    operations_research::MPVariable *var = solver.MakeIntVar(0.0, 1.0, col_name);
			       variables.push_back(var);
			    var->SetInteger(true);                                   // on définit la variable comme binaire
			    solver.MutableObjective()->SetCoefficient(var, 0.0);     // définir le coefficient de l'objectif pour cette variable à 0
		    }
	    }

	    //==============================================================================================================================================================================================================================================================================================
	   /* // color unicity per pixel/contraintes d'unicité de la couleur par pixel
	    for (int p = 0; p < nbVert; p++) {                                 // parcours chaque pixel
		    std::string row_name = "colorunicity_" + std::to_string(p);     // construit un nom unoqie  pour chaque contrainte de couleur unie par pixel
		    glp_set_row_name(lp, irow, row_name.c_str());                   // assigne le nom de la ligne irow dans le pb glpk, avec irow l'index de la ligne
		    glp_set_row_bnds(lp, irow, GLP_FX, 1, 0);     // définit les bornes(limites) de la ligne irow dans le bp glpk. la contrainte est fixée à 1 donc le somme des variables de couleur pour ce pixel égale à 1

		    for (int imat = 0; imat < nbMat; imat++) {     // parcours chaque matériau
			    ia[innz] = irow;                            // définit l'index de la ligne irow pour la position innz dans le tableau ia
			    ja[innz] = p * nbMat + imat;     // définit l'index de la colonne pour la position innz dans le tableau ja. on calcul l'index de la colonne en utilisant : p*nbMat+imat, ce qui représente la variable de couleur pour le pixel p et le matériau imat
			    ar[innz] = 1.;                   // définit la valeur du coeff de la matrice de contraintes à 1 pour la position innz dans le tableau ar

			    innz++;     // on incrémente le compteur innz pour la prochaine position dans les tableau ia,ja et ar
		    }

		    irow++;     // on incrémente également le compteur de lignes irow pour la prochaine contrainte de couleur unie
	    }
	    */
	    // code convertit
	    // contraintes de couleur d'unicité par pixel

	    for (int p = 0; p < nbVert; p++) {
		    std::string row_name = "colorunicity_" + std::to_string(p);
		    operations_research::MPConstraint *const constraint = solver.MakeRowConstraint(1.0, 1.0, row_name);
		    for (int imat = 0; imat < nbMat; imat++) {
			    constraint->SetCoefficient(variables[p * nbMat + imat], 1.0);
		    }
	    }

	    //=================================================================================================================================================================================================================================================================================================

	    // vf per coarse mixed cell per material/contraintes du pourcentage de la fraction volumique par cellule
	    const int nbCells_source = ACellsIDs_source->getNbElems();
	    /*// on obtient le nb total de cellules sources
	    for (int i = 0; i < nbCells_source; i++) {                                 // parcours chaque cellule source.
		    kmds::TCellID cid = ACellsIDs_source->get(i);                           // cid obtinent l'id de la cellule source
		    if ((*AVarMixedCells_source)[cid]) {                                    // on verifie si la cellule cid est une cellule mixte
			    kmds::TCellID subCid_first = (*AVarOldCells2firstSubCells)[cid];     // si oui , subCid_first obtient l'id de la 1ère sous-cellule associée à cid

			    int sum_vf2int = 0;                            // on initialise une somme pour vérifier les fractions volumiques
			    for (int imat = 0; imat < nbMat; imat++) {     // parcours chaque matériau imat
				    std::string row_name =
				       "vfpreserve_" + std::to_string(cid) + "_"
				       + std::to_string(
				          imat);     // définit le nom de la ligne pour la contrainte de préservation volumique et l'associe à 'index de la ligne irow dans GLPK
				    glp_set_row_name(lp, irow, row_name.c_str());

				    int vf2int = ANbSubPixels * Afp_source->getFracPres(imat, cid);     // calcul la fraction volumique pour le matériau imat dans la cellule cid
				    glp_set_row_bnds(lp, irow, GLP_FX, vf2int,
				                     0);     // définit une contrainte d'égalité (GLP_FIX) pour que la somme des pixels assignés à ce matériau soit égale à vf2int

				    for (int isub = 0; isub < ANbSubPixels; isub++) {                       // boucle for pour chaque sous-pixels isub de la cellule mixte.
					    kmds::TCellID vert = (*AVarCells2vertices)[subCid_first + isub];     // vert obtient l’identifiant du sommet associé au sous-pixel.

					    ia[innz] = irow;     // Ajoute les coefficients pour les contraintes dans les tableaux ia, ja et ar.
					    ja[innz] = vert * nbMat + imat;
					    ar[innz] = 1.;

					    innz++;     // on incrémente le compteur innz pour le prochain coeff
				    }

				    sum_vf2int += vf2int;     // on incrémente sum_vf2int de la valeur volumique courante
				    irow++;                   // on incrémente également le compteur de lignes irow pour la prochaine contrainte
			    }

			    // check correctness for vf * nbPixels treated as an integer
			    if (sum_vf2int != ANbSubPixels) {
				    std::cout << "sum_vf2int " << sum_vf2int << std::endl;
				    throw kmds::KException(
				       "InterfaceNodesPosSmoothVF_buildSubGraph_glpk_2D : sum_vf2int not equal to number of pixels.");     // Vérifie que sum_vf2int est égal au nombre de sous-pixels (ANbSubPixels). Si ce n’est pas le cas, lance une exception.
			    }
		    }
	    }
*/
	    // code convertit
	    //  =======================================
	    //  vf per coarse mixed cell per material/contraintes du pourcentage de la fraction volumique par cellule

	    for (int i = 0; i < nbCells_source; i++) {                                 // parcours chaque cellule source.
		    kmds::TCellID cid = ACellsIDs_source->get(i);                           // cid obtinent l'id de la cellule source
		    if ((*AVarMixedCells_source)[cid]) {                                    // on verifie si la cellule cid est une cellule mixte
			    kmds::TCellID subCid_first = (*AVarOldCells2firstSubCells)[cid];     // si oui , subCid_first obtient l'id de la 1ère sous-cellule associée à cid
			    int sum_vf2int = 0;                                                  // on initialise une somme pour vérifier les fractions volumiques
			    for (int imat = 0; imat < nbMat; imat++) {                           // parcours chaque matériau imat
				    std::string row_name = "vfpreserve_" + std::to_string(cid) + "_" + std::to_string(imat);
				    operations_research::MPConstraint *C0 =
				       solver.MakeRowConstraint(ANbSubPixels * Afp_source->getFracPres(imat, cid), ANbSubPixels * Afp_source->getFracPres(imat, cid), row_name);
				    int vf2int = ANbSubPixels * Afp_source->getFracPres(imat, cid);
				    for (int isub = 0; isub < ANbSubPixels; isub++) {
					    kmds::TCellID vert = (*AVarCells2vertices)[subCid_first + isub];
					    C0->SetCoefficient(variables[vert * nbMat + imat], 1.0);
				    }
				    sum_vf2int += vf2int;
			    }
			    if (sum_vf2int != ANbSubPixels) {
				    std::cout << "sum_vf2int " << sum_vf2int << std::endl;
				    throw kmds::KException("InterfaceNodesPosSmoothVF_buildSubGraph_glpk_2D : sum_vf2int not equal to number of pixels.");
			    }
		    }
	    }

	    //======================================================================================================================================================================================================================================================================================================================================
/*
	    // neighboring similarity/ définition des contraintes pour la Similarité des couleurs voisines
	    for (int p = 0; p < nbVert; p++) {     // Boucle sur chaque pixel p et chaque matériau imat
		    for (int imat = 0; imat < nbMat; imat++) {
			    std::string col_pos_name =
			       "neighbourcolor_pos_" + std::to_string(p) + "_"
			       + std::to_string(imat);     // Définit les noms et les bornes des colonnes pour les variables de similarité de couleur positive
			    glp_set_col_name(lp, icol, col_pos_name.c_str());
			    glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);

			    glp_set_obj_coef(lp, icol, 1.);     // Définit les coefficients des fonctions objectifs pour ces variables.

			    icol++;     // on incrémente également le compteur de colonnes icol pour la prochaine contrainte

			    std::string col_neg_name =
			       "neighbourcolor_neg_" + std::to_string(p) + "_"
			       + std::to_string(imat);     // Définit les noms et les bornes des colonnes pour les variables de similarité de couleur négative
			    glp_set_col_name(lp, icol, col_neg_name.c_str());
			    glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);

			    glp_set_obj_coef(lp, icol, 1.);     // Définit les coefficients des fonctions objectifs pour ces variables.

			    icol++;     // on incrémente également le compteur de colonnes icol pour la prochaine contrainte

			    std::string row_name =
			       "neighbourcolor_abs_" + std::to_string(p) + "_"
			       + std::to_string(imat);     // Définit le nom de la ligne et les bornes pour la contrainte de similarité absolue des couleurs.
			    glp_set_row_name(lp, irow, row_name.c_str());
			    glp_set_row_bnds(lp, irow, GLP_FX, 0, 0);

			    // Cpk/Ajoute les coefficients pour les contraintes Cpk.
			    ia[innz] = irow;
			    ja[innz] = p * nbMat + imat;
			    ar[innz] = 1.;
			    innz++;

			    // De/Ajoute les coefficients pour les contraintes  De .
			    ia[innz] = irow;
			    ja[innz] = nbVert * nbMat + 2 * (p * nbMat + imat) + 0;
			    ar[innz] = 1.;
			    innz++;     // Incrémente innz pour chaque coefficient ajouté.

			    // de/Ajoute les coefficients pour les contraintes des tableaux ia, ja et ar.
			    ia[innz] = irow;
			    ja[innz] = nbVert * nbMat + 2 * (p * nbMat + imat) + 1;
			    ar[innz] = -1.;
			    innz++;     // Incrémente innz pour chaque coefficient ajouté.

			    Kokkos::View<kmds::TCellID *> neighbors;
			    int nb = AGraph->getNeighbors(p, neighbors);     // Obtient les voisins du pixel p.

			    for (int in = 0; in < neighbors.size(); in++) {     // Boucle pour ajouter les coefficients des voisins dans les tableaux ia, ja et ar.

				    ia[innz] = irow;
				    ja[innz] = neighbors[in] * nbMat + imat;
				    ar[innz] = -(1. / neighbors.size());

				    innz++;     // Incrémente innz pour la prochaine contrainte.
			    }

			    irow++;     // on incrémente irow pour la prochaine contrainte
		    }
	    }
	    //============================================================================================================================================================================================
	    */for (int p = 0; p < nbVert; p++) {
		    for (int imat = 0; imat < nbMat; imat++) {
			    std::string col_pos_name = "neighbourcolor_pos_" + std::to_string(p) + "_" + std::to_string(imat);
			    operations_research::MPVariable *var_pos = solver.MakeNumVar(0.0, solver.infinity(), col_pos_name);
			    solver.MutableObjective()->SetCoefficient(var_pos, 1.0);
			    std::string col_neg_name = "neighbourcolor_neg_" + std::to_string(p) + "_" + std::to_string(imat);
			    operations_research::MPVariable *var_neg = solver.MakeNumVar(0.0, solver.infinity(), col_neg_name);
			    solver.MutableObjective()->SetCoefficient(var_neg, 1.0);
			    std::string row_name = "neighbourcolor_abs_" + std::to_string(p) + "_" + std::to_string(imat);
			    operations_research::MPConstraint *C1 = solver.MakeRowConstraint(0.0, 0.0, row_name);
			    C1->SetCoefficient(variables[p * nbMat + imat], 1.0);
			    C1->SetCoefficient(variables[nbVert * nbMat + 2 * (p * nbMat + imat) + 0], 1.0);
			    C1->SetCoefficient(variables[nbVert * nbMat + 2 * (p * nbMat + imat) + 1], -1.0);
			    Kokkos::View<kmds::TCellID *> neighbors;
			    int nb = AGraph->getNeighbors(p, neighbors);     // Obtient les voisins du pixel p.

			    for (int in = 0; in < neighbors.size(); in++) {
				    C1->SetCoefficient(variables[neighbors[in] * nbMat + imat], -(1.0 / neighbors.size()));
			    }
		    }

	    }
	    //========================================================================================================================================================================================================================================================
/*
	    // known color for subpixels of pure cells/ définition des contraintes des pixels purs
	    for (int p = 0; p < nbVert; p++) {     // Boucle sur chaque pixel p pour vérifier s’il est dans vert2puremat.
		    if (vert2puremat.find(p) != vert2puremat.end()) {
			    int imat = vert2puremat[p];     // Définit les contraintes pour les pixels purs en associant les pixels aux matériaux connus.
			    std::string row_name = "purepixels_" + std::to_string(p);
			    glp_set_row_name(lp, irow, row_name.c_str());
			    glp_set_row_bnds(lp, irow, GLP_FX, 1, 0);

			    ia[innz] = irow;
			    ja[innz] = p * nbMat + imat;
			    ar[innz] = 1.;

			    innz++;

			    irow++;
		    }
	    }
*/
	    for (int p = 0; p < nbVert; p++) {     // Boucle sur chaque pixel p pour vérifier s’il est dans vert2puremat.
		    if (vert2puremat.find(p) != vert2puremat.end()) {
			    int imat = vert2puremat[p];     // Définit les contraintes pour les pixels purs en associant les pixels aux matériaux connus.
			    std::string row_name = "purepixels_" + std::to_string(p);
			    operations_research::MPConstraint *C2 = solver.MakeRowConstraint(1.0, 1.0, row_name);
			    C2->SetCoefficient(variables[p * nbMat + imat], 1.0);
		    }

	    }

	   /* //===========================================================================================================================================================================================================================================================================
	    // Affiche les informations sur les lignes, colonnes et coefficients non nuls.
	    std::cout << "row " << nbRows << " irow " << irow << std::endl;
	    std::cout << "col " << nbCols << " icol " << icol << std::endl;
	    std::cout << "nnz " << nnz << " innz " << innz << std::endl;

	    // add one because that is the indices starting index
	    // Incrémente les indices des colonnes car GLPK utilise des indices de 1.
	    for (int j = 1; j <= nnz; j++) {
		    ja[j]++;
	    }
*/
	    // convertit
	    std::cout <<"row"<< solver.NumConstraints() <<" irow " << irow << std::endl;
	    std::cout << "col" << solver.NumVariables() << " icol " << icol << std::endl;
	    std::cout << "nnz " << nnz << " innz " << innz << std::endl;

	    operations_research::MPSolver::ResultStatus result_status = solver.Solve();
	    // Check that the problem has an optimal solution.
	    if (result_status != operations_research::MPSolver::OPTIMAL) {
		    std::cerr << "le problème n'a pas de solution optimal." << std::endl;
	    }
	    std::cout << "Objective value =" << solver.Objective().Value() << std::endl;



/*

//================================================================================================================================================================================================================================================================================================
//Charge la matrice de contraintes dans GLPK.
        glp_load_matrix(lp, nnz, ia, ja, ar);

        // glp_term_out( GLP_OFF );
	     //active la sortie terminale
	     //initialise les paramètres pour GLPK
	     //Active la présolve
        glp_term_out(GLP_ON);//active la sortie terminale

        glp_iocp glpParams;//on déclare la structure nommée glpParams
        glp_init_iocp(&glpParams);//  initialise la structure glpParams
        glpParams.presolve = GLP_ON;//Active la présolve

        glpParams.tm_lim = 300000;//on definit la limite de temps à 300 000 millisecondes dans glpParams
//        glpParams.ps_heur = GLP_ON;

        glp_write_lp(lp, NULL, "cplex.txt"); //on ecrit le problème sous format LP de CPLEX

        glp_write_mps(lp, GLP_MPS_FILE, NULL, "mps.txt"); //on ecrit le problème sous format MPS de CPLEX ce qui nous permet le debogage ou l'examen du pb de programmation linéaire

	     int glpErr = 0; // glperr capture le code d'erreur du solveur
        glpErr = glp_intopt(lp, &glpParams); // glp_intop résout le pb en utilisant la programmation linéire en nombres entiers avec les paramètres spécifiés
        switch (glpErr) { //verifie le résultat de l'optimisation, si glpErr est 0 alors l'optimisation est réussi.Sinon ,il y a eu un pb lors de la résolution
            case 0:
                std::cout << "GLP OK" << std::endl;
                break;
            default:
                std::cout << "pb solving in GLP." << std::endl;
                break;
        }

        glpErr = glp_mip_status(lp); // vérifie le status de la résolution en utilisant glp_mip_status. le statut peut indiquer si la solution est non définis(GLP_UNDEF),si la solution est optimal(GLP_OPT) et si elle est faisable (GLP_FEAS).si le statut est indéfini ou s'il n'y a pas de solution faisable, une exception est lancée
        switch (glpErr) {
            case GLP_UNDEF:
                std::cout << " MIP solution is undefined" << std::endl;
                throw kmds::KException("SubMapping::boundaryDiscretization MIP solution is undefined");
                break;
            case GLP_OPT:
                std::cout << " MIP solution is integer optimal" << std::endl;
                break;
            case GLP_FEAS:
                std::cout << " MIP solution is integer feasible" << std::endl;
                break;
            case GLP_NOFEAS:
                std::cout << "problem has no integer feasible solution" << std::endl;
                throw kmds::KException("SubMapping::boundaryDiscretization problem has no integer feasible solution");
                break;
            default:
                throw kmds::KException("SubMapping::boundaryDiscretization glp_intopt unknown return code.");
        }

        glp_print_mip(lp, "pixelAssignment.txt");// on imprime la solution de l'optimisation dans un fichier pixelAssignment.txt

        for (int i = 0; i < nbCells_target; i++) { //boucle pour chaque cellule cible
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid]; // pour chaque cellule, récupère le sommet(vert) associé
            if (vert != kmds::NullID) {// si le sommet est valide (cad vert != NullID), on détermine le matériau(mat) en vérifiant les valeurs des colonnes du modèle GLPK
                int mat = nbMat; // default;
                bool found = false;
                for(int imat=0; imat<nbMat; imat++) {
                    int color = glp_mip_col_val(lp, vert * nbMat + imat + 1);
                    if(color == 1) { // si on trouve une couleur pour un matériau alors cette couleur est assignée à la cellule
                        if(found) {
                            std::cout<<"QQQQQQQQQQQQQQQQ"<<std::endl;
                        }
                        mat = imat;
                        found = true;
                    }
                }
                Ama_target->setMaterial(mat, cid); // mise à jour du matériau de la cellule cible avec Ama_target-> setMaterial(mat,cid)
            } else {
//                Ama_target->setMaterial(nbMat, cid);
            }
        }
*/

	     //converit
	     operations_research::MPSolverParameters param;
	     param.SetIntegerParam(operations_research::MPSolverParameters::PRESOLVE,operations_research::MPSolverParameters::PRESOLVE_ON);
	     solver.SetTimeLimit(absl::Duration 300000);
	     if (result_status !=operations_research::MPSolver::OPTIMAL && result_status!= operations_research::MPSolver::FEASIBLE){
		     std::cout << "Le probleme n'a pas de solution optimale." << std::endl;
	     } // a mettre en commentaires puis tester code en haut
	     switch (result_status) {
	     case operations_research::MPSolver::OPTIMAL:
		     std::cout << "MIP solution est un entier" << std::endl;
		     break;
	     case operations_research::MPSolver:: FEASIBLE:
		     std::cout<<"La solution MIP  est un entier, est faisable "<< std::endl;
		     break;
	     default;
		      throw kmds::KException("SubMapping::boundaryDiscretization glp_intopt unknown return code.");
	     }
	     std::ofstream solution_file("pixelAssignment.txt");
	     if(!solution_file){
		     std::cout<<"No such file";
	     }
	     if(solution_file.is_open()){
		     for(int i=0;i<solver.NumVariables();i++){
			     solution_file<<"variable"<<i<<":Value="<<solver.variable(i)-> solution_value()<<std::endl;
		     }
		     solution_file.close();
	     }

	     for (int i = 0; i < nbCells_target; i++) { //boucle pour chaque cellule cible
		     kmds::TCellID cid = ACellsIDs_target->get(i);

		     kmds::TCellID vert = (*AVarCells2vertices)[cid]; // pour chaque cellule, récupère le sommet(vert) associé
		     if (vert != kmds::NullID) {// si le sommet est valide (cad vert != NullID), on détermine le matériau(mat) en vérifiant les valeurs des colonnes du modèle GLPK
			     int mat = nbMat; // default;
			     bool found = false;
			     for(int imat=0; imat<nbMat; imat++) {
				     int color = solver.variable(vert*nbMat + imat)->solution_value();
				     if(color == 1) { // si on trouve une couleur pour un matériau alors cette couleur est assignée à la cellule
					     if(found) {
						     std::cout<<"QQQQQQQQQQQQQQQQ"<<std::endl;
					     }
					     mat = imat;
					     found = true;
				     }
			     }
			     Ama_target->setMaterial(mat, cid); // mise à jour du matériau de la cellule cible avec Ama_target-> setMaterial(mat,cid)
		     } else {
			     //                Ama_target->setMaterial(nbMat, cid);
		     }
	     }





	     /*
//=================================================================================================================================================================================================================================================================================================================
//libère la mémoire allouée pour les tableau ia,ja et ar
        free(ia);
        free(ja);
        free(ar);
//supprime le pb GLPK POUR LIBÉRÉR LES RESSOURCES Associées
        glp_delete_prob(lp);

    }

    /*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
