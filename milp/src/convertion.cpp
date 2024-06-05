//
// Created by guirassye on 17/05/24.
//

//#include "nbVert.h"
#include "ortools/linear_solver/linear_solver.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace operations_research;
 from ortools.linear_solver import pywraplp ;// on importe le module pywraplp d'Or-Tools qui sera utilisé pour la création et la résolution de pb de programmation linéaire et mixte
   nbVert=AGraph.getNbVec();// on récupère le nb de sommet dans le graphe AGraph
 nbMat= Afp_source.getNbMaterials();//récupère nb de matériel l'objet Afp_source
 nbCells_target=ACellsIDs_target.getNbElems();//on récupère le nombre de cellules cibles dans l'objet ACellsIDs_target
 vert2puremat= {};//on va créer une map qui associe chaque sommet(kmds::TCellID) à un entier représentant le matériau pur.
 for i in range(nbCells_target):cid = ACellsIDs_target.get(i) // boucle qui parcours chaque cellule cible
         vert = AVarCells2vertices[cid];//récupère l'id de la cellule cible actuelle
   if vert != kmds.NullID:// vérifie si le sommet est valide
     cid_old=AVarSubCells2oldCells[cid]// récupère l'ancienne cellule associée
   if not Afp_source.isMixedCells(cid_old):// vérifie si cette ancienne cellule n'est pas mixte
         for imat in range(nbMat)://parcourt chaque matériau
             if Afp_source.getFracPres(imat,cid_old)==1.0: //vérifie si le matériau imat a une fraction de cellule cid_old. cela signifie que imat est le matériau pur de cette cellule
                vert2puremat[vert]=imat; // si la condition prcdt est remplie, ajoute l'association entre le sommet vert et le matériau imat à la map vert2puremat
             break;
         int nbMixedCells_source = 0; // initialise  un compteur pour le nombre de cellules mixtes
         for i in range(ACellsIDs_source.getNbElems()): // utilise une réduction en parallèle pour compter les cellules / récupère le nombre de cellules sources
              cid=ACellsIDs_source.get(i): // récupère l'id de la cellule source actuelle
              if AVarmixedCells_source[cid]:// incrémente nbMixedCells_source si la cellule est mixte
                 nbMixedCells_source += 1;


         // =======================================
         // ORTOOLS
         // =======================================
         //Initialisation de GLPK/on cherche aussi à minimiser la fonction objective
         def solve_optimization_problem(nbVert,nbMat,nbMixedCells_source,ANbSubpixels,vert2puremat,AGraph): //définir une fonction en puthon pour créer le pb d'optimisation
           solver= pywraplp.Solver.CreateSolver('SCIP'); // Céer un solveur SCIP avec OR-TOOLS
                    if not solver: // vérifier si le solveur est disponible, sinon, afficher un message et retourner
                          print("SCIP solver not availabe.")
                    return;
         solver.Objective().SetMinimization(); // définir l'objectif du solveur comme une minimisation
         //Détermination du nb total de variables et de contraintes
         nbRows=0;//initialise le compteur du nb de lignes à 0
         nbCols=0;// initialise le compteur de nb de colonnes à 0
         nnz = 0;// initialise le compteur du nb d'éléments non nuls dans la matrice des contraintes à 0

         // pixel color variables/variables de couleur des pixels
         nbCols += nbVert * nbMat; //ajoute au compteur de colonnes le nb tot de variables de couleur des pixels, soit nbVert*nbMat

         // color unicity per pixel/unicité de la couleur par pixel: 1 pixel ne peut avoir qu'une seule couleur attribuée
         nbRows += nbVert;// ajout au compteur  de ligne une contrainte pour chaque pixel, -> chaque pixel n'a qu'une seule couleure
         nnz += nbVert * nbMat; // ajoute au compteur d'éléments != 0 le produit du nb de pixels par le nb de matériaux

         // vf per coarse mixed cell	/ calcul du pourcentage de matériaux  par cellule mixte: la somme des fractions volumiques des matériaux dans une seule cellule est égale à un nombre spécifique.
         nbRows += nbMixedCells_source * nbMat; // ajoute au compteur lignes des contraintes pour les fractions volumiques de chaque cellule mixte et chaque matériaux
         nnz += nbMixedCells_source * nbMat *ANbSubPixels;

         // neighboring similarity/ cette contrainte assure que les couleurs des pixels voisins sont les mêmes
         nbRows += nbVert * nbMat; // ajoute au compteur de ligne des contraintes pour la similarité de voisinge entre les pixels.
         nbCols += 2 * nbVert * nbMat; // ajoute au compteur des colonnes des variables supplémentaires pour les contraintes de voisinage
         nnz += 3 * nbVert * nbMat + nbMat * AGraph.getNbEdges(); // ajoute au compteur d'éléments des contraintes de similarité de voisinage et des arêtes du graphe des pixels

         // ajouts des contraintes des variables dans GLPK
         // pure pixels
         nbRows += vert2puremat.size(); // ajoute au compteur de lignes un contrinte pour chaque pixel pur( où la couleur est fixée à un matériau spécifique)
         nnz += vert2puremat.size(); // ajoute au compteur d'éléments non nuls les contributions de chaque pixel pur.

         // pixel color variables/contraintes de couleur
         // represented as booleans
         //créer un dictionnaire pour stocker les variables
         pixelColor_Variables ={};

         for p in range(nbVert):
           for imat in range(nbMat):
               col_name=f"pixelColor_{p}_{imat}" ; //on crée un nom unique pour chaque variable de couleur de pixel
         pixelColor_Variables[col_name] = solver.BoolVar(col_name); // crée une variable binaire avec le nom donné (col_name)
         solver.Objective().setCoefficient(pixelColor_Variables[col_name],0.0); // on définie le coefficient de l'objectif pour cette variable à 0
         for p in range(nbVert): // parcours chaque pixel
           solver.Add(sum(pixelColor_Variable[f"pixelColor_{p}_{imat}"]);// on ajoute une contrainte pour chaque pixel pour s'assurer qu'il ne peut avoir qu'une seule couleur
                          for imat in range(nbMat))==1):
                              status = solver.Solve; // résout le pb d'optimisation et stocke le statut de la solution
         if status == pywraplp.Solver.OPTIMAL: print("Solution trouvée") // vérifie si une solution optimale a été trouvée
                                      for p in range(nbVert): // parcours chaque pixel
                                          for imat in range(nbMat):// parcours chaque matériau
                                              varName=f"pixelColor_{p}_{imat}"; // on nomme la variable
                                              print(f'{varName} = {pixelColor_Variables[varName].solution_value()}'); //affiche la valeur de solution pour chaque variable de ouleur de pixel
         else:
              print("Aucune solution optimale trouvée");

         // color unicity per pixel/contraintes d'unicité de la couleur par pixel
         for p in range(nbVert): //  parcours chaque pixel
           row_name=f"colorUnicity_{p}"; // construit un nom pour chaque contrainte de couleur unie par pixel
         constraint = solver.RowConstraint(1,1,row_name); // définit la contrainte avec une somme égale à 1
         for imat in range(nbMat): // parcours chaque matériau
           varName=f"pixelColor_{p}_{imat}" ; // on nomme la variable
         constraint.SetCoefficient(pixelColor_Variables[varName],1); // définit le coefficient de la matrice de contrainte à 1
         status = solver.Solve(); // résolution du problème
// affichage des résultats
         if status == pywraplp.Solver.OPTIMAL: print("Solution trouvée") // vérifie si une solution optimale a été trouvée
                                      for p in range(nbVert): // parcours chaque pixel
                                          for imat in range(nbMat):// parcours chaque matériau
                                              varName=f"pixelColor_{p}_{imat}"; // on nomme la variable
         print(f'{varName} = {pixelColor_Variables[varName].solution_value()}'); //affiche la valeur de solution pour chaque variable de couleur de pixel
         else:
           print("Aucune solution optimale trouvée");

         // vf per coarse mixed cell per material/contraintes du pourcentage de la fraction volumique par cellule
         nbCells_Source = ACellsIDs_source.getNbElems(); // on obtient le nb total de cellules sources
         for i in range(nbCells_source): //parcours chaque cellule source.
           cid = ACellsIDs_source.get(i); //cid obtinent l'id de la cellule source
           if AVarMixedCells_source[cid]:
             subCid_first = AVar01dCells2firstSubCells[cid];
             sum_vf2int = 0;
             for imat in range(nbMat):
               row_name = f"vfpreserve"_{cid}_{imat};
               vf2int = ANbSubPixels * Afp_source.getFracPres(imat,cid);
               constraint = solver.RowConstraint(vf2int,vf2int,row_name);
               for isub in range(ANbSubPixels):
                 vert = AVarCells2vertices[subCid_first + isub]
                 varName= f"pixelColor"_{vert}_{imat};
                 constraint.SetCoefficient(pixelColor_Variables[varName],1);
                 sum_vf2int +=vf2int;
             if sum_vf2int != ANbSubPixels:
               print(f"sum_vf2int{sum_vf2int}");
               raise Exception("sum_vf2int n'est pas égale au nombre de pixels voulu."); // Signale une erreur en interrompant l'exécution du programme et en affichant un message d'erreur

//=================================================================================================================================================================================================================================================================================================================================

               //code convertit en C++
               // =======================================
               // ORTOOLS
               // =======================================
               //Initialisation de GLPK/on cherche aussi à minimiser la fonction objective
               operations_research::MPSolver solver("InterfaceNodesPosSmoothVF", operations_research::MPSolver::SCIP_MIXED_INTEGER_PROGRAMMING);
               if( ! solver) {
	               std::cerr << "echec dans la création du solver." << std::endl;
	               return;
               }
               // compute problem size
               //Détermination du nb total de variables et de contraintes
               int nbRows = 0; //initialise le compteur du nb de lignes à 0
               int nbCols = 0; // initialise le compteur de nb de colonnes à 0
               int nnz = 0; // initialise le compteur du nb d'éléments non nuls dans la matrice des contraintes à 0

               // pixel color variables/variables de couleur des pixels
               nbCols += nbVert * nbMat; //ajoute au compteur de colonnes le nb tot de variables de couleur des pixels, soit nbVert*nbMat

               // color unicity per pixel/unicité de la couleur par pixel: 1 pixel ne peut avoir qu'une seule couleur attribuée
               nbRows += nbVert; // ajout au compteur  de ligne une contrainte pour chaque pixel, -> chaque pixel n'a qu'une seule couleure
               nnz += nbVert * nbMat; // ajoute au compteur d'éléments != 0 le produit du nb de pixels par le nb de matériaux

               // vf per coarse mixed cell	/ calcul du pourcentage de matériaux  par cellule mixte: la somme des fractions volumiques des matériaux dans une seule cellule est égale à un nombre spécifique.
               nbRows += nbMixedCells_source * nbMat; // ajoute au compteur lignes des contraintes pour les fractions volumiques de chaque cellule mixte et chaque matériaux
               nnz += nbMixedCells_source * nbMat *ANbSubPixels; // ajoute au compteur d'éléments !=0 les contributions de chaque cellule mixte et de chaque matériau multiplié par le nombre de sius-pixels

               // neighboring similarity/ cette contrainte assure que les couleurs des pixels voisins sont les mêmes
               nbRows += nbVert * nbMat; // ajoute au compteur de ligne des contraintes pour la similarité de voisinge entre les pixels.
               nbCols += 2 * nbVert * nbMat; // ajoute au compteur des colonnes des variables supplémentaires pour les contraintes de voisinage
               nnz += 3 * nbVert * nbMat + nbMat * AGraph->getNbEdges(); // ajoute au compteur d'éléments des contraintes de similarité de voisinage et des arêtes du graphe des pixels

               // ajouts des contraintes des variables dans GLPK
               // pure pixels
               nbRows += vert2puremat.size(); // ajoute au compteur de lignes un contrinte pour chaque pixel pur( où la couleur est fixée à un matériau spécifique)
               nnz += vert2puremat.size(); // ajoute au compteur d'éléments non nuls les contributions de chaque pixel pur.

               // pixel color variables/contraintes de couleur
               // represented as booleans
               for(int p = 0;p<nbVert;p++){
	               for(int imat=0;imat<nbMat;imat++){
		               std::string col_name = "pixelColor_" + std::to_string(p)+"_"+std::to_string(imat);
		               operations_research::MPVariable* var = solver.MakeIntVar(0.0,1.0,col_name);
		               variables.push_back(var);
		               var-> SetInteger(true);// on définit la variable comme binaire
		               solver.MutableObjective()->SetCoefficient(var,0.0); //définir le coefficient de l'objectif pour cette variable à 0
	               }
               }

               // contraintes de couleur d'unicité par pixel
               for(int p = 0;p<nbVert;p++){
	               std::string row_name = "colorunicity_" + std::to_string(p);
	               operations_research::MPConstraint*const constraint = solver.MakeRowConstraint(1.0,1.0,row_name);
	               for(int imat=0;imat<nbMat;imat++){
		               constraint->SetCoefficient(variables[p*nbMat+imat],1.0);
	               }
               }

               // vf per coarse mixed cell per material/contraintes du pourcentage de la fraction volumique par cellule
               const int nbCells_source = ACellsIDs_source->getNbElems(); // on obtient le nb total de cellules sources
               for (int i = 0; i < nbCells_source; i++) {                 // parcours chaque cellule source.
	               kmds::TCellID cid = ACellsIDs_source->get(i);           // cid obtinent l'id de la cellule source
	               if ((*AVarMixedCells_source)[cid]) {                    // on verifie si la cellule cid est une cellule mixte
		               kmds::TCellID subCid_first =
		                  (*AVarOldCells2firstSubCells)[cid];         // si oui , subCid_first obtient l'id de la 1ère sous-cellule associée à cid
		               int sum_vf2int = 0;                            // on initialise une somme pour vérifier les fractions volumiques
		               for (int imat = 0; imat < nbMat; imat++) {     // parcours chaque matériau imat
			               std::string row_name = "vfpreserve_" + std::to_string(cid) + "_" + std::to_string(imat);
			               operations_research::MPConstraint *c = solver.MakeRowConstraint(ANbSubPixels * Afp_source->getFracPres(imat, cid),ANbSubPixels * Afp_source->getFracPres(imat, cid), row_name);
			               int vf2int = ANbSubPixels * Afp_source->getFracPres(imat, cid);
			               for (int isub = 0; isub < ANbSubPixels; isub++) {
				               kmds::TCellID vert = (*AVarCells2vertices)[subCid_first + isub];
				               c->SetCoefficient(variables[vert * nbMat + imat], 1.0);
			               }
			               sum_vf2int += vf2int;
		               }
		               if (sum_vf2int != ANbSubPixels) {
			               std::cout << "sum_vf2int " << sum_vf2int << std::endl;
			               throw kmds::KException("InterfaceNodesPosSmoothVF_buildSubGraph_glpk_2D : sum_vf2int not equal to number of pixels.");
		               }
	               }
               }

               // neighboring similarity/ définition des contraintes pour la Similarité des couleurs voisines
               for(int p=0; p<nbVert; p++){
	               for(int imat=0;imat<nbMat;imat++){
		               std::string col_pos_name = "neighbourcolor_pos_" + std::to_string(p) + "_" + std::to_string(imat);
		               operations_research::MPVariable* var_pos = solver.MakeNumVar(0.0,solver.infinity(), col_pos_name);
		               solver.MutableObjective()->SetCoefficient(var_pos,1.0);
		               std::string col_neg_name= "neighbourcolor_neg_" + std::to_string(p) + "_" + std::to_string(imat);
		               operations_research::MPVariable* var_neg = solver.MakeNumVar(0.0,solver.infinity(), col_neg_name);
		               solver.MutableObjective()->SetCoefficient(var_neg,1.0);
		               std::string row_name = "neighbourcolor_abs_" + std::to_string(p) + "_" + std::to_string(imat);
		               operations_research::MPConstraint* C1 = solver.MakeRowConstraint(0.0,0.0, row_name);
		               C1->SetCoefficient(variables[p*nbMat +imat],1.0);
		               C1->SetCoefficient(variables[nbVert*nbMat+2*(p*nbMat+imat)+0],1.0);
		               C1->SetCoefficient(variables[nbVert*nbMat+2*(p*nbMat+imat)+1],-1.0);
		               Kokkos::View<kmds::TCellID *> neighbors;
		               int nb = AGraph->getNeighbors(p, neighbors); //Obtient les voisins du pixel p.

		               for (int in=0; in<neighbors.size(); in++) {
			               C1-> SetCoefficient(variables[neighbors[in] *nbMat + imat],-(1.0 / neighbors.size()));
		               }

	               }
	               solver.Solve();
               }

               // known color for subpixels of pure cells/ définition des contraintes des pixels purs
               for(int p=0; p<nbVert; p++) {     // Boucle sur chaque pixel p pour vérifier s’il est dans vert2puremat.
	               if (vert2puremat.find(p) != vert2puremat.end()) {
		               int imat = vert2puremat[p];     // Définit les contraintes pour les pixels purs en associant les pixels aux matériaux connus.
		               std::string row_name = "purepixels_" + std::to_string(p);
		               operations_research::MPConstraint* C2 = solver.MakeRowConstraint(1.0,1.0, row_name);
		               C2->SetCoefficient(variables[p*nbMat +imat],1.0);
	               }
	               solver.Solve();
               }



               // Affiche les informations sur les lignes, colonnes et coefficients non nuls.
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
//MPSolver solver("InterfaceNodesPosSmoothVF_assignMIP_xD",MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
if(!solver) {
	LOG(FATAL) << "SCIP solver unavailable.";
}

//Conversion des variables GLPK en OR-TOOLS
const int nbVert = AGraph->getNBVec(); //nb de sommets
const int nbMat = Afp_source->getNBMaterials(); //Nombre de
//Déclaration des variables
std::vector<std::vector<MPVariable*>>pixel_color(nbVert,std::vector<MPVariable*>(nbMat));
for (int i = 0;i < nbVert;i++){
	for(int imat= 0; imat < nbMat; imat++)
	{
		pixel_color[i][imat] = solver.MakeBoolVar("pixelcolor_"+std::to_string(i) + "_" + std ::to_string(imat));
	}
}

//Contrainte : unicité de couleur par pixel
for(int i = 0;i < nbVert;i++){
	LinearExpr sum;
	for(int imat = 0; imat<nbMat;imat++){
		sum += pixel_color[i][imat];
	}
	solver.MakeRowConstraint(sum == 1.0,"colorunicity_",+ std::to_string(i));
}
// contrainte




