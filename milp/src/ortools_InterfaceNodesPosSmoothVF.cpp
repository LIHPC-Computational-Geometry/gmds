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
                                           kmds::Graph *AGraph,
                                           const int ANbSubPixels)
{
	// =======================================
	// Initialisation des variables, avec nbVert= nb de sommets dans un graph , nbMat = nb de matériaux, nbCells_target= nb de cellules cibles
	const int nbVert = AGraph->getNbVec();              // récupère nb de sommets dans le graph Agraph
	const int nbMat = Afp_source->getNbMaterials();     // récupère nb de matériel l'objet Afp_source

	const int nbCells_target = ACellsIDs_target->getNbElems();     // on récupère le nombre de cellules cibles dans l'objet ACellsIDs_target

	// ici on va créer une map pour les pixels purs.Le but est d'associer chaque pixel à un matériau pur. Grâce à cette étape on peut savoir quels pixels sont purs et à quels matériaux ils sont associés

	std::map<kmds::TCellID, int> vert2puremat;     // on va créer une map qui associe chaque sommet(kmds::TCellID) à un entier représentant le matériau pur.
	for (int i = 0; i < nbCells_target; i++) {     // boucle qui parcours chaque cellule cible
		kmds::TCellID cid = ACellsIDs_target->get(i);     // récupère l'id de la cellule cible actuelle

		kmds::TCellID vert = (*AVarCells2vertices)[cid];              // récupère le sommet associé à cette cellule
		if (vert != kmds::NullID) {                                   // vérifie si le sommet est valide
			kmds::TCellID cid_old = (*AVarSubCells2oldCells)[cid];     // récupère l'ancienne cellule associée
			if (!Afp_source->isMixedCell(cid_old)) {                   // vérifie si cette ancienne cellule n'est pas mixte
				for (int imat = 0; imat < nbMat; imat++) {              // parcourt chaque matériau
					if (Afp_source->getFracPres(imat, cid_old)
					    == 1.) {     // vérifie si le matériau imat a une fraction de cellule cid_old. cela signifie que imat est le matériau pur de cette cellule
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
	   ACellsIDs_source->getNbElems(),            // utilise une réduction en parallèle pour compter les cellules / récupère le nombre de cellules sources
	   KOKKOS_LAMBDA(const int i, int &sum) {     // la fonction lambda qui prend en indice i et une réfécrence à sum
		   kmds::TCellID cid = ACellsIDs_source->get(i);     // récupère l'id de la cellule source actuelle
		   if ((*AVarMixedCells_source)[cid]) {              // incrémente sum si la cellule est mixte
			   sum++;
		   }
	   },
	   nbMixedCells_source);
	std::cout << "Initialisation des données effectuée" << std::endl;

	// code convertit
	//  =======================================
	//  ORTOOLS
	//  =======================================
	// Initialisation de GLPK/on cherche aussi à minimiser la fonction objective

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

	int irow = 1;     // start at 1, not 0 !
	int icol = 1;     // start at 1, not 0 !
	int innz = 1;     // start at 1, not 0 !

	// pixel color variables/variables de couleur des pixels
	nbCols += nbVert * nbMat;     // ajoute au compteur de colonnes le nb tot de variables de couleur des pixels, soit nbVert*nbMat

	// color unicity per pixel/unicité de la couleur par pixel: 1 pixel ne peut avoir qu'une seule couleur attribuée
	nbRows += nbVert;          // ajout au compteur  de ligne une contrainte pour chaque pixel, -> chaque pixel n'a qu'une seule couleure
	nnz += nbVert * nbMat;     // ajoute au compteur d'éléments != 0 le produit du nb de pixels par le nb de matériaux

	// vf per coarse mixed cell	/ calcul du pourcentage de matériaux  par cellule mixte: la somme des fractions volumiques des matériaux dans une seule cellule est égale à un nombre spécifique.
	nbRows +=
	   nbMixedCells_source * nbMat;     // ajoute au compteur lignes des contraintes pour les fractions volumiques de chaque cellule mixte et chaque matériaux
	nnz += nbMixedCells_source * nbMat * ANbSubPixels;     // ajoute au compteur d'éléments !=0 les contributions de chaque cellule mixte et de chaque matériau multiplié par le nombre de sius-pixels

	// neighboring similarity/ cette contrainte assure que les couleurs des pixels voisins sont les mêmes
	nbRows += nbVert * nbMat;         // ajoute au compteur de ligne des contraintes pour la similarité de voisinge entre les pixels.
	nbCols += 2 * nbVert * nbMat;     // ajoute au compteur des colonnes des variables supplémentaires pour les contraintes de voisinage
	nnz += 3 * nbVert * nbMat
	       + nbMat * AGraph->getNbEdges(0);     // ajoute au compteur d'éléments des contraintes de similarité de voisinage et des arêtes du graphe des pixels

	// ajouts des contraintes des variables dans GLPK
	// pure pixels
	nbRows += vert2puremat.size();     // ajoute au compteur de lignes un contrinte pour chaque pixel pur( où la couleur est fixée à un matériau spécifique)
	nnz += vert2puremat.size();        // ajoute au compteur d'éléments non nuls les contributions de chaque pixel pur.

	// code convertit
	//  =======================================
	//  ORTOOLS
	//  =======================================
	std::vector<operations_research::MPVariable *> variables;
	operations_research::MPSolver solver("InterfaceNodesPosSmoothVF", operations_research::MPSolver::SCIP_MIXED_INTEGER_PROGRAMMING);
	for (int p = 0; p < nbVert; p++) {     // boule parcourue pour chaque pixel et matériau
		for (int imat = 0; imat < nbMat; imat++) {
			std::string col_name =
			   "pixelColor_" + std::to_string(p) + "_" + std::to_string(imat);                // crée un nom unique pour chaque variable de couleur bde pixel
			operations_research::MPVariable *var = solver.MakeIntVar(0.0, 1.0, col_name);     // définit le nom de la variable dans le problème GLPK
			variables.push_back(var);     // ajoute un élément à la fin du vecteur variables ce qui augmente la taille du vecteur de un à chaque ajout(le nouvel élément est placé à la fin du vecteur)
			var->SetInteger(true);                                   // on définit la variable comme binaire
			solver.MutableObjective()->SetCoefficient(var, 0.0);     // définir le coefficient de l'objectif pour cette variable à 0
		}
	}
	std::cout << "Variables ajoutées avec succès" << std::endl;
	for (int p = 0; p < nbVert; p++) {                                 // parcours chaque pixel
		std::string row_name = "colorunicity_" + std::to_string(p);     // construit un nom unique  pour chaque contrainte de couleur unie par pixel
		operations_research::MPConstraint *const constraint =
		   solver.MakeRowConstraint(1.0, 1.0, row_name);     // création de la contrainte linéaire  et difnit les bornes
		if (constraint == nullptr) {                         // on vérifie si la contrainte a bien été créée
			std::cerr << "Failed to create constraint for" << row_name << std::endl;
			return;
		}
		for (int imat = 0; imat < nbMat; imat++) {     // parcours chaque matériau
			constraint->SetCoefficient(variables[p * nbMat + imat], 1.0);
		}
	}
	std::cout << "Contraintes d'unicité ajoutées avec succès" << std::endl;

	//=================================================================================================================================================================================================================================================================================================

	// vf per coarse mixed cell per material/contraintes du pourcentage de la fraction volumique par cellule
	const int nbCells_source = ACellsIDs_source->getNbElems();
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
				if (C0 == nullptr) {     // on vérifie si la contrainte a bien été créée
					std::cerr << "Failed to create constraint for" << row_name << std::endl;
					return;
				}
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
	std::cout << "Contraintes du pourcentage de la fraction volumique par cellule ajoutées avec succès" << std::endl;

	//============================================================================================================================================================================================

	operations_research::MPObjective *const objective = solver.MutableObjective();
	for (int p = 0; p < nbVert; p++) {
		for (int imat = 0; imat < nbMat; imat++) {
			std::string col_pos_name = "neighbourcolor_pos_" + std::to_string(p) + "_" + std::to_string(imat);
			operations_research::MPVariable *var_pos = solver.MakeNumVar(0.0, solver.infinity(), col_pos_name);
			variables.push_back(var_pos);
			solver.MutableObjective()->SetCoefficient(var_pos, 1.0);
			std::string col_neg_name = "neighbourcolor_neg_" + std::to_string(p) + "_" + std::to_string(imat);
			operations_research::MPVariable *var_neg = solver.MakeNumVar(0.0, solver.infinity(), col_neg_name);
			variables.push_back(var_neg);
			solver.MutableObjective()->SetCoefficient(var_neg, 1.0);
			std::string row_name = "neighbourcolor_abs_" + std::to_string(p) + "_" + std::to_string(imat);
			operations_research::MPConstraint *C1 = solver.MakeRowConstraint(0.0, 0.0, row_name);
			if (C1 == nullptr) {     // on vérifie si la contrainte a bien été créée
				std::cerr << "Failed to create constraint for" << row_name << std::endl;
				return;
			}
			C1->SetCoefficient(variables[p * nbMat + imat], 1.0);
			C1->SetCoefficient(variables[nbVert * nbMat + 2 * (p * nbMat + imat) + 0], 1.0);
			C1->SetCoefficient(variables[nbVert * nbMat + 2 * (p * nbMat + imat) + 1], -1.0);
			Kokkos::View<kmds::TCellID *> neighbors;
			int nb = AGraph->getNeighbors(p, neighbors);

			for (int in = 0; in < neighbors.size(); in++) {
				C1->SetCoefficient(variables[neighbors[in] * nbMat + imat], -(1.0 / neighbors.size()));
			}
		}
	}

	objective->SetMinimization();
	std::cout << "Contraintes d'ajout des voisins ajoutées avec succès" << std::endl;
	for (int p = 0; p < nbVert; p++) {     // Boucle sur chaque pixel p pour vérifier s’il est dans vert2puremat.
		if (vert2puremat.find(p) != vert2puremat.end()) {
			int imat = vert2puremat[p];     // Définit les contraintes pour les pixels purs en associant les pixels aux matériaux connus.
			std::string row_name = "purepixels_" + std::to_string(p);
			operations_research::MPConstraint *C2 = solver.MakeRowConstraint(1.0, 1.0, row_name);
			C2->SetCoefficient(variables[p * nbMat + imat], 1.0);
		}
	}
	std::cout << "Contraintes de vérifications des pixels purs  ajoutées avec succès" << std::endl;

	// convertit
	std::cout << "row" << solver.NumConstraints() << " irow " << irow << std::endl;
	std::cout << "col" << solver.NumVariables() << " icol " << icol << std::endl;
	std::cout << "nnz " << nnz << " innz " << innz << std::endl;

	solver.SetTimeLimit(absl::Milliseconds(600000));
	operations_research::MPSolver::ResultStatus result_status = solver.Solve();
	// Check that the problem has an optimal solution.
	if (result_status != operations_research::MPSolver::OPTIMAL) {
		std::cerr << "le problème n'a pas de solution optimal." << std::endl;
	}
	std::cout << "Objective value =" << solver.Objective().Value() << std::endl;

	// converit
	operations_research::MPSolverParameters param;
	param.SetIntegerParam(operations_research::MPSolverParameters::PRESOLVE, operations_research::MPSolverParameters::PRESOLVE_ON);
	if (result_status != operations_research::MPSolver::OPTIMAL && result_status != operations_research::MPSolver::FEASIBLE) {
		std::cout << "Le probleme n'a pas de solution optimale." << std::endl;
	}     // a mettre en commentaires puis tester code en haut

	switch (result_status) {
	case operations_research::MPSolver::OPTIMAL: std::cout << "MIP solution est un entier" << std::endl; break;
	case operations_research::MPSolver::FEASIBLE: std::cout << "La solution MIP est un entier, est faisable" << std::endl; break;
	case operations_research::MPSolver::INFEASIBLE: std::cerr << "Le problème est infaisable." << std::endl; break;
	case operations_research::MPSolver::NOT_SOLVED: std::cerr << "Le problème n'a pas été résolu." << std::endl; break;
	default: throw std::runtime_error("SubMapping::boundaryDiscretization unknown return code.");
	}

	if (result_status == operations_research::MPSolver::OPTIMAL || result_status == operations_research::MPSolver::FEASIBLE) {
		// Continuez avec votre logique actuelle.
		std::ofstream solution_file("pixelAssignment.txt");
		if (!solution_file) {
			std::cout << "No such file";
		}
		if (solution_file.is_open()) {
			for (int i = 0; i < solver.NumVariables(); i++) {
				solution_file << "variable" << i << ":Value=" << solver.variable(i)->solution_value() << std::endl;
			}
			solution_file.close();
		}

		for (int i = 0; i < nbCells_target; i++) {
			kmds::TCellID cid = ACellsIDs_target->get(i);

			kmds::TCellID vert = (*AVarCells2vertices)[cid];
			if (vert != kmds::NullID) {
				int mat = nbMat;     // default
				bool found = false;
				for (int imat = 0; imat < nbMat; imat++) {
					int color = solver.variable(vert * nbMat + imat)->solution_value();
					if (color == 1) {
						if (found) {
							std::cout << "QQQQQQQQQQQQQQQQ" << std::endl;
						}
						mat = imat;
						found = true;
					}
				}
				Ama_target->setMaterial(mat, cid);
			}
		}
		std::cout << "Assignation des matériaux effectuée" << std::endl;
	}
	else {
		// Gérer les cas où le problème n'est pas résolu de manière optimale ou faisable.
		std::cerr << "Le problème n'a pas été résolu de manière optimale ou faisable." << std::endl;
	}
}
}