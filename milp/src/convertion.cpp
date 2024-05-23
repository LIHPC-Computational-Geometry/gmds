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

// =======================================
// OR-TOOLS
// =======================================
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



// =======================================
// GLPK code explication
// =======================================

void
InterfaceNodesPosSmoothVF_assignMIP_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
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

	//Initialisation des variables, avec nbVert= nb de sommets dans un graph , nbMat = nb de matériaux, nbCells_target= nb de cellules cibles
	const int nbVert = AGraph->getNbVec(); // récupère nb de sommets dans le graph Agraph
	const int nbMat = Afp_source->getNbMaterials(); //récupère nb de matériel l'objet Afp_source


	const int nbCells_target = ACellsIDs_target->getNbElems();// récupère le nb cellules cibles dans l'objet

	//ici on va créer une map pour les pixels purs.Le but est d'associer chaque pixel a un matériau pur. Grâce à cette étape on peut savoir quels pixels sont purs et à quels matériau ils sont associés
	std::map<kmds::TCellID, int> vert2puremat;
	for (int i = 0; i < nbCells_target; i++) {
		kmds::TCellID cid = ACellsIDs_target->get(i);

		kmds::TCellID vert = (*AVarCells2vertices)[cid];
		if (vert != kmds::NullID) {
			kmds::TCellID cid_old = (*AVarSubCells2oldCells)[cid];
			if(!Afp_source->isMixedCell(cid_old)) {
				for(int imat=0; imat<nbMat; imat++) {
					if(Afp_source->getFracPres(imat, cid_old) == 1.) {
						vert2puremat.emplace(vert, imat);
						break;
					}
				}
			}
		}
	}

	// calcul du nombre de cellules mixtes dans le graphe. cellules mixtes = cellules qui possèdeent plrs matériaux
	int nbMixedCells_source = 0;
	Kokkos::parallel_reduce(ACellsIDs_source->getNbElems(),
	                        KOKKOS_LAMBDA(const int i, int& sum) {
		                        kmds::TCellID cid = ACellsIDs_source->get(i);
		                        if ((*AVarMixedCells_source)[cid]) {
			                        sum++;
		                        }
	                        },
	                        nbMixedCells_source);


	// =======================================
	// GLPK
	// =======================================
// Initialisation de GLPK/on cherche aussi à minimiser la fonction objective
	glp_prob *lp;
	int *ia, *ja;
	double *ar;

	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MIN);

	// compute problem size
	//Détermination du nb total de variables et de contraintes
	int nbRows = 0;
	int nbCols = 0;
	int nnz = 0;

	// pixel color variables /variables de couleur des pixels
	nbCols += nbVert * nbMat;

	// color unicity per pixel/unicité de la couleur par pixel: 1 pixel ne peut avoir qu'une seule couleur attribuée
	nbRows += nbVert;
	nnz += nbVert * nbMat;
	Format_5x5_2D.txt 1
	   // vf per coarse mixed cell/ calcul du pourcentage de matériaux  par cellule mixte: la somme des fractions volumiques des matériaux dans une seule cellule est égale à un nombre spécifique.
	   nbRows += nbMixedCells_source * nbMat;
	nnz += nbMixedCells_source * nbMat *ANbSubPixels;

	// neighboring similarity/ cette contrainte assure que les couleurs des pixels voisins sont les mêmes
	nbRows += nbVert * nbMat;
	nbCols += 2 * nbVert * nbMat;
	nnz += 3 * nbVert * nbMat + nbMat * AGraph->getNbEdges();

	// pure pixels
	nbRows += vert2puremat.size();
	nnz += vert2puremat.size();

// ajouts des contraintes des variables dans GLPK
	// enforce constraints
	glp_add_rows(lp, nbRows);
	glp_add_cols(lp, nbCols);

	ia = (int *) calloc(1 + nnz, sizeof(int));
	ja = (int *) calloc(1 + nnz, sizeof(int));
	ar = (double *) calloc(1 + nnz, sizeof(double));


	int irow = 1; // start at 1, not 0 !
	int icol = 1; // start at 1, not 0 !
	int innz = 1; // start at 1, not 0 !

	// pixel color variables/contraintes de couleur
	// represented as booleans
	for(int p=0; p<nbVert; p++) {
		for(int imat=0; imat<nbMat; imat++) {
			std::string col_name = "pixelcolor_" + std::to_string(p) + "_" + std::to_string(imat);
			glp_set_col_name(lp, icol, col_name.c_str());
			glp_set_col_kind(lp, icol, GLP_BV);
			// glp_set_col_bnds(lp, icol, GLP_DB, 0, 1); // not necessary because set by GLP_BV

			glp_set_obj_coef(lp, icol, 0.);

			icol++;
		}
	}

	// color unicity per pixel/contraintes d'unicité de la couleur par pixel
	for(int p=0; p<nbVert; p++) {
		std::string row_name = "colorunicity_" + std::to_string(p);
		glp_set_row_name(lp, irow, row_name.c_str());
		glp_set_row_bnds(lp, irow, GLP_FX, 1, 0);

		for(int imat=0; imat<nbMat; imat++) {
			ia[innz] = irow;
			ja[innz] = p * nbMat + imat;
			ar[innz] = 1.;

			innz++;
		}

		irow++;
	}

	// vf per coarse mixed cell per material/contraintes du pourcentage de la fraction volumique par cellule
	const int nbCells_source = ACellsIDs_source->getNbElems();
	for (int i = 0; i < nbCells_source; i++) {
		kmds::TCellID cid = ACellsIDs_source->get(i);
		if ((*AVarMixedCells_source)[cid]) {
			kmds::TCellID subCid_first = (*AVarOldCells2firstSubCells)[cid];

			int sum_vf2int = 0;
			for(int imat=0; imat<nbMat; imat++) {
				std::string row_name = "vfpreserve_" + std::to_string(cid) + "_" + std::to_string(imat);
				glp_set_row_name(lp, irow, row_name.c_str());

				int vf2int = ANbSubPixels * Afp_source->getFracPres(imat, cid);
				glp_set_row_bnds(lp, irow, GLP_FX, vf2int, 0);

				for (int isub = 0;
				     isub < ANbSubPixels; isub++) {
					kmds::TCellID vert = (*AVarCells2vertices)[subCid_first + isub];

					ia[innz] = irow;
					ja[innz] = vert * nbMat + imat;
					ar[innz] = 1.;

					innz++;
				}

				sum_vf2int += vf2int;
				irow++;
			}

			// check correctness for vf * nbPixels treated as an integer
			if(sum_vf2int != ANbSubPixels) {
				std::cout<<"sum_vf2int "<<sum_vf2int<<std::endl;
				throw kmds::KException("InterfaceNodesPosSmoothVF_buildSubGraph_glpk_2D : sum_vf2int not equal to number of pixels.");
			}
		}
	}

	// neighboring similarity / définition des contraintes pour la Similarité des couleurs voisines
	for(int p=0; p<nbVert; p++) {
		for(int imat=0; imat<nbMat; imat++) {
			std::string col_pos_name = "neighbourcolor_pos_" + std::to_string(p) + "_" + std::to_string(imat);
			glp_set_col_name(lp, icol, col_pos_name.c_str());
			glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);

			glp_set_obj_coef(lp, icol, 1.);

			icol++;

			std::string col_neg_name = "neighbourcolor_neg_" + std::to_string(p) + "_" + std::to_string(imat);
			glp_set_col_name(lp, icol, col_neg_name.c_str());
			glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);

			glp_set_obj_coef(lp, icol, 1.);

			icol++;

			std::string row_name = "neighbourcolor_abs_" + std::to_string(p) + "_" + std::to_string(imat);
			glp_set_row_name(lp, irow, row_name.c_str());
			glp_set_row_bnds(lp, irow, GLP_FX, 0, 0);

			// Cpk
			ia[innz] = irow;
			ja[innz] = p * nbMat + imat;
			ar[innz] = 1.;
			innz++;

			// De
			ia[innz] = irow;
			ja[innz] = nbVert * nbMat + 2 * (p * nbMat + imat) + 0;
			ar[innz] = 1.;
			innz++;

			// de
			ia[innz] = irow;
			ja[innz] = nbVert * nbMat + 2 * (p * nbMat + imat) + 1;
			ar[innz] = - 1.;
			innz++;

			Kokkos::View<kmds::TCellID *> neighbors;
			int nb = AGraph->getNeighbors(p, neighbors);

			for (int in=0; in<neighbors.size(); in++) {

				ia[innz] = irow;
				ja[innz] = neighbors[in] * nbMat + imat;
				ar[innz] = - (1. / neighbors.size());

				innz++;
			}

			irow++;
		}
	}

	// known color for subpixels of pure cells/ définition des contraintes des pixels purs
	for(int p=0; p<nbVert; p++) {
		if(vert2puremat.find(p) != vert2puremat.end()) {
			int imat = vert2puremat[p];
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


	std::cout<<"row "<<nbRows<<" irow "<<irow<<std::endl;
	std::cout<<"col "<<nbCols<<" icol "<<icol<<std::endl;
	std::cout<<"nnz "<<nnz<<" innz "<<innz<<std::endl;

	// add one because that is the indices starting index
	for(int j=1; j<=nnz; j++) {
		ja[j]++;
	}

	glp_load_matrix(lp, nnz, ia, ja, ar);

	// glp_term_out( GLP_OFF );
	glp_term_out(GLP_ON);

	glp_iocp glpParams;
	glp_init_iocp(&glpParams);
	glpParams.presolve = GLP_ON;

	glpParams.tm_lim = 300000;
	//        glpParams.ps_heur = GLP_ON;

	glp_write_lp(lp, NULL, "cplex.txt");

	glp_write_mps(lp, GLP_MPS_FILE, NULL, "mps.txt");

	int glpErr = 0;
	glpErr = glp_intopt(lp, &glpParams);
	switch (glpErr) {
	case 0:
		std::cout << "GLP OK" << std::endl;
		break;
	default:
		std::cout << "pb solving in GLP." << std::endl;
		break;
	}

	glpErr = glp_mip_status(lp);
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

	glp_print_mip(lp, "pixelAssignment.txt");

	for (int i = 0; i < nbCells_target; i++) {
		kmds::TCellID cid = ACellsIDs_target->get(i);

		kmds::TCellID vert = (*AVarCells2vertices)[cid];
		if (vert != kmds::NullID) {
			int mat = nbMat; // default;
			bool found = false;
			for(int imat=0; imat<nbMat; imat++) {
				int color = glp_mip_col_val(lp, vert * nbMat + imat + 1);
				if(color == 1) {
					if(found) {
						std::cout<<"QQQQQQQQQQQQQQQQ"<<std::endl;
					}
					mat = imat;
					found = true;
				}
			}
			Ama_target->setMaterial(mat, cid);
		} else {
			//                Ama_target->setMaterial(nbMat, cid);
		}
	}

	free(ia);
	free(ja);
	free(ar);

	glp_delete_prob(lp);

}
 */