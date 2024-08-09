/*----------------------------------------------------------------------------*/
//
// Created by calderans on 29/06/2023.
//
/*----------------------------------------------------------------------------*/
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include <fstream>
#include "gmds/claire/CGNSWriter3D.h"
#include <iomanip>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
using namespace claire;

#define Structured CG_Structured
#define RealDouble CG_RealDouble
#define BCTypeUserDefined CG_BCTypeUserDefined
#define PointRange CG_PointRange
#define Unstructured CG_Unstructured
#define TRI_3 CG_TRI_3
/*----------------------------------------------------------------------------*/
CGNSWriter3D::CGNSWriter3D(Blocking3D *ABlocking)
  :m_blocks(ABlocking),m_mesh(nullptr)
{}
/*----------------------------------------------------------------------------*/
CGNSWriter3D::CGNSWriter3D(Mesh *AMesh)
  :m_mesh(AMesh),m_blocks(nullptr)
{}
/*----------------------------------------------------------------------------*/
CGNSWriter3D::CGNSWriter3D()
  :m_mesh(nullptr)
{
	m_blocks = new Mesh(MeshModel(DIM3 | N | F | R | N2F | N2R | F2N | F2R | R2N | R2F));
}
/*----------------------------------------------------------------------------*/
CGNSWriter3D::~CGNSWriter3D()
{}
/*----------------------------------------------------------------------------*/
void CGNSWriter3D::initialize(const std::string &AOutFileName, const std::string &dir){

	std::cout<<"========================================================="<<std::endl;
	std::cout<<"Initialize CGNS Writer"<<std::endl;
	std::cout<<"========================================================="<<std::endl;

	m_stream = new std::ofstream(AOutFileName, std::ios::out);
	if (!m_stream){
		std::string s ="Impossible to create a CGNS File: "+AOutFileName;
		throw GMDSException(s);
	}

	int n = AOutFileName.length() + dir.length();

	std::string filename;
	filename.resize(n+1);
	filename = dir + AOutFileName;

	if(cg_open(filename.c_str(),CG_MODE_WRITE,&m_indexFile) != CG_OK){
		std::cout<<"Section error : "<<cg_get_error()<<std::endl;
	}

	char basename[AOutFileName.length()];
	strcpy(basename, AOutFileName.c_str());

	m_cellDim = 3;
	m_physdim = 3;

	cg_base_write(m_indexFile,basename,m_cellDim,m_physdim,&m_indexBase);

}
/*----------------------------------------------------------------------------*/
void CGNSWriter3D::writeZones()
{

	std::cout << "=========================================================" << std::endl;
	std::cout << "Start writing zones" << std::endl;
	std::cout << "=========================================================" << std::endl;

	int index_tf = 0;
	int id_bc = 0;

	Variable<int>* farfield = m_blocks->getVariable<int,GMDS_NODE>("Farfield");
	Variable<int>* paroi = m_blocks->getVariable<int,GMDS_NODE>("Paroi");
	Variable<int>* out = m_blocks->getOrCreateVariable<int,GMDS_NODE>("Sortie");
	Variable<int>* solid = m_blocks->getOrCreateVariable<int,GMDS_REGION>("Solide");


	for (auto bid : m_blocks->regions()) {
		Region b = m_blocks->get<Region>(bid);

		//std::cout << "Writing zone " << b.id()+1 << std::endl;

		// Le nom de la zone correspond au nom du bloc, dans notre cas on l'appel juste
		//  "Block_<Block_ID>"
		std::stringstream ss;
		ss << (solid->value(b.id()) == 1 ? "SOLIDE" : "FLUIDE ") << std::setw(5) << std::setfill('0') << b.id()+1;     // Pour le moment on utilise les faces du mesh comme bloc
		std::string name = ss.str();
		char zonename[32];
		strcpy(zonename, name.c_str());

		int discrI = VarDiscrI->value(bid);
		int discrJ = VarDiscrJ->value(bid);
		int discrK = VarDiscrK->value(bid);

		// Ici on va remplir les infos sur le nombre de sommets, mailles etc du bloc
		cgsize_t zone_size[9] = {discrI, discrJ, discrK, discrI - 1, discrJ - 1, discrK - 1,0,0,0};

		cg_zone_write(m_indexFile, m_indexBase, zonename, zone_size, Structured, &m_indexZone);

		int size = discrI * discrJ * discrK;
		auto x_coords = new double[size];
		auto y_coords = new double[size];
		auto z_coords = new double[size];

		std::vector<TCellID> b_grid = m_block_grid->value(b.id());

		int i_coord = 0;
		for (int i = 0; i < discrI*discrJ*discrK; i++) {
			Node node = m_blocks->get<Node>(b_grid[i]);
			x_coords[i_coord] = node.X();
			y_coords[i_coord] = node.Y();
			z_coords[i_coord] = node.Z();
			i_coord++;
		}

		int index_coord;

		char coord_nameX[32];
		strcpy(coord_nameX, "CoordinateX");
		char coord_nameY[32];
		strcpy(coord_nameY, "CoordinateY");
		char coord_nameZ[32];
		strcpy(coord_nameZ, "CoordinateZ");

		cg_coord_write(m_indexFile, m_indexBase, m_indexZone, RealDouble, coord_nameX, x_coords, &index_coord);
		cg_coord_write(m_indexFile, m_indexBase, m_indexZone, RealDouble, coord_nameY, y_coords, &index_coord);
		cg_coord_write(m_indexFile, m_indexBase, m_indexZone, RealDouble, coord_nameZ, z_coords, &index_coord);


		//Writing Boundary Condition
		for(int iFace = 0; iFace<6; iFace++){
			Face face = b.get<Face>()[iFace];

			if(face.get<Region>().size() == 2){
				int id_voisin = face.getIDs<Region>()[0] == b.id() ? face.getIDs<Region>()[1] : face.getIDs<Region>()[0];

				char interface[33];
				std::stringstream interface_ss;
				interface_ss << "S" << std::setw(5) << std::setfill('0') << face.id()+1 << " (F" << std::setw(5) << std::setfill('0') << b.id()+1 << " & F" << std::setw(5)
				             << std::setfill('0') << id_voisin+1 << ")";
				std::string interface_s = interface_ss.str();
				strcpy(interface, interface_s.c_str());
				//std::cout << "\t -> connection " << interface_s << std::endl;

				std::stringstream name2_ss;
				name2_ss << "FLUIDE " << std::setw(5) << std::setfill('0') << id_voisin+1;
				std::string name2 = name2_ss.str();     // Pour le moment on utilise les faces du mesh comme bloc
				char zonename2[32];
				strcpy(zonename2, name2.c_str());

				cgsize_t pts1[9], pts2[9];
				int idmin1,idinter1,idmax1;

				int discrI_vois = VarDiscrI->value(id_voisin),discrJ_vois = VarDiscrJ->value(id_voisin),discrK_vois = VarDiscrK->value(id_voisin);

				switch (iFace) {
				case 0:
					pts1[0] = 1;	// X1
					pts1[1] = 1;	// Y1
					pts1[2] = 1;	// Z1

					pts1[3] = discrI;	// X2
					pts1[4] = 1;		// Y2
					pts1[5] = 1;		// Z2

					pts1[6] = discrI;	// X3
					pts1[7] = discrJ;	// Y3
					pts1[8] = 1;		// Z3

					idmin1 = 0;
					idinter1 = 1;
					idmax1 = 2;
					break;
				case 2:
					pts1[0] = 1;
					pts1[1] = 1;
					pts1[2] = 1;

					pts1[3] = discrI;
					pts1[4] = 1;
					pts1[5] = 1;

					pts1[6] = discrI;
					pts1[7] = 1;
					pts1[8] = discrK;

					idmin1 = 0;
					idinter1 = 1;
					idmax1 = 5;
					break;
				case 5:
					pts1[0] = discrI;
					pts1[1] = 1;
					pts1[2] = 1;

					pts1[3] = discrI;
					pts1[4] = discrJ;
					pts1[5] = 1;

					pts1[6] = discrI;
					pts1[7] = discrJ;
					pts1[8] = discrK;

					idmin1 = 1;
					idinter1 = 2;
					idmax1 = 6;
					break;
				case 3:
					pts1[0] = 1;
					pts1[1] = discrJ;
					pts1[2] = 1;

					pts1[3] = discrI;
					pts1[4] = discrJ;
					pts1[5] = 1;

					pts1[6] = discrI;
					pts1[7] = discrJ;
					pts1[8] = discrK;

					idmin1 = 3;
					idinter1 = 2;
					idmax1 = 6;
					break;
				case 4:
					pts1[0] = 1;
					pts1[1] = 1;
					pts1[2] = 1;

					pts1[3] = 1;
					pts1[4] = discrJ;
					pts1[5] = 1;

					pts1[6] = 1;
					pts1[7] = discrJ;
					pts1[8] = discrK;

					idmin1 = 0;
					idinter1 = 3;
					idmax1 = 7;
					break;
				case 1:
					pts1[0] = 1;	// X1
					pts1[1] = 1;	// Y1
					pts1[2] = discrK;	// Z1

					pts1[3] = discrI;	// X2
					pts1[4] = 1;		// Y2
					pts1[5] = discrK;	// Z2

					pts1[6] = discrI;	// X3
					pts1[7] = discrJ;	// Y3
					pts1[8] = discrK;	// Z3

					idmin1 = 4;
					idinter1 = 5;
					idmax1 = 6;
					break;
				default: break;
				}

				int id_glob_n0 = b.getIDs<Node>()[idmin1];
				int id_glob_n1 = b.getIDs<Node>()[idmax1];
				int id_glob_inter1 = b.getIDs<Node>()[idinter1];

				Region voisin = m_blocks->get<Region>(id_voisin);

				std::vector<int> indices2;
				indices2.resize(3);

				std::vector<TCellID> nid_vois = voisin.getIDs<Node>();
				for(int i = 0; i<8; i++){
					if(nid_vois[i] == id_glob_n0){
						indices2[0] = i;
					}else if(nid_vois[i] == id_glob_inter1){
						indices2[1] = i;
					}else if(nid_vois[i] == id_glob_n1){
						indices2[2] = i;
					}
				}


				for(int i = 0; i<3; i++) {
					switch (indices2[i]) {
					case 0:
						pts2[0+(3*i)] = 1;     // X1
						pts2[1+(3*i)] = 1;     // Y1
						pts2[2+(3*i)] = 1;     // Z1
						break;
					case 1:
						pts2[0+(3*i)] = discrI_vois;
						pts2[1+(3*i)] = 1;
						pts2[2+(3*i)] = 1;
						break;
					case 2:
						pts2[0+(3*i)] = discrI_vois;
						pts2[1+(3*i)] = discrJ_vois;
						pts2[2+(3*i)] = 1;
						break;
					case 3:
						pts2[0+(3*i)] = 1;
						pts2[1+(3*i)] = discrJ_vois;
						pts2[2+(3*i)] = 1;
						break;
					case 4:
						pts2[0+(3*i)] = 1;
						pts2[1+(3*i)] = 1;
						pts2[2+(3*i)] = discrK_vois;
						break;
					case 5:
						pts2[0+(3*i)] = discrI_vois;
						pts2[1+(3*i)] = 1;
						pts2[2+(3*i)] = discrK_vois;
						break;
					case 6:
						pts2[0+(3*i)] = discrI_vois;
						pts2[1+(3*i)] = discrJ_vois;
						pts2[2+(3*i)] = discrK_vois;
						break;
					case 7:
						pts2[0+(3*i)] = 1;
						pts2[1+(3*i)] = discrJ_vois;
						pts2[2+(3*i)] = discrK_vois;
						break;
					default: break;
					}
				}
				int transform[3];


				bool filtre[3] = {false, false, false};
				bool filtre_vois[3] = {false, false, false};

				int isize1[3] = {discrI,discrJ,discrK};
				int isize2[3] = {discrI_vois,discrJ_vois,discrK_vois};

				int ind1, ind2;
				int val1, val2;
				int ipnts1[3],ipnts2[3];
				for (int k=0; k<3; k++){
					ipnts1[k] = pts1[k];
					ipnts2[k] = pts1[6+k];
				}
				_getIndicesIdAndVal(ipnts1,ipnts2,filtre,ind1,val1);
				filtre[ind1] = true;

				for (int k=0; k<3; k++){
					ipnts1[k] = pts2[k];
					ipnts2[k] = pts2[6+k];
				}
				_getIndicesIdAndVal(ipnts1,ipnts2,filtre_vois,ind2,val2);
				filtre_vois[ind2] = true;

				bool is_val1_min;
				if (val1 == 1)
					is_val1_min = true;
				else if (isize1[ind1] == val1)
					is_val1_min = false;
				else {
					std::string	message = "Erreur dans cg_1to1_write step 1,  val1 = "+ std::to_string(val1) + " n'est ni au min ni au max " + std::to_string(isize1[ind1]);
					throw GMDSException (message);
				}

				// idem avec val2
				bool is_val2_min;
				if (val2 == 1)
					is_val2_min = true;
				else if (isize2[ind2] == val2)
					is_val2_min = false;
				else {
					std::string	message = "Erreur dans cg_1to1_write step 1,  val2 = "+ std::to_string(val2) + " n'est ni au min ni au max " + std::to_string(isize2[ind2]);
					throw GMDSException (message);
				}

				transform[ind1] = ind2+1;
				if (is_val1_min == is_val2_min)
					transform[ind1] = -transform[ind1];




				for (int k=0; k<3; k++){
					ipnts1[k] = pts1[k];
					ipnts2[k] = pts1[3+k];
				}
				_getIndicesIdAndVal(ipnts1,ipnts2,filtre,ind1,val1);
				filtre[ind1] = true;

				for (int k=0; k<3; k++){
					ipnts1[k] = pts2[k];
					ipnts2[k] = pts2[3+k];
				}
				_getIndicesIdAndVal(ipnts1,ipnts2,filtre_vois,ind2,val2);
				filtre_vois[ind2] = true;


				if (val1 == 1)
					is_val1_min = true;
				else if (isize1[ind1] == val1)
					is_val1_min = false;
				else {
					std::string	message = "Erreur dans cg_1to1_write step 2,  val1 = "+ std::to_string(val1) + " n'est ni au min ni au max " + std::to_string(isize1[ind1]);
					throw GMDSException (message);
				}

				// idem avec val2
				if (val2 == 1)
					is_val2_min = true;
				else if (isize2[ind2] == val2)
					is_val2_min = false;
				else {
					std::string	message = "Erreur dans cg_1to1_write step 2,  val2 = "+ std::to_string(val2) + " n'est ni au min ni au max " + std::to_string(isize2[ind2]);
					throw GMDSException (message);
				}

				transform[ind1] = ind2+1;
				if (is_val1_min != is_val2_min)
					transform[ind1] = -transform[ind1];

				for (int k=0; k<3; k++){
					ipnts1[k] = pts1[3+k];
					ipnts2[k] = pts1[6+k];
				}
				_getIndicesIdAndVal(ipnts1,ipnts2,filtre,ind1,val1);
				filtre[ind1] = true;

				for (int k=0; k<3; k++){
					ipnts1[k] = pts2[3+k];
					ipnts2[k] = pts2[6+k];
				}
				_getIndicesIdAndVal(ipnts1,ipnts2,filtre_vois,ind2,val2);
				filtre_vois[ind2] = true;


				if (val1 == 1)
					is_val1_min = true;
				else if (isize1[ind1] == val1)
					is_val1_min = false;
				else {
					std::string	message = "Erreur dans cg_1to1_write step 3,  val1 = "+ std::to_string(val1) + " n'est ni au min ni au max " + std::to_string(isize1[ind1]);
					throw GMDSException (message);
				}

				// idem avec val2
				if (val2 == 1)
					is_val2_min = true;
				else if (isize2[ind2] == val2)
					is_val2_min = false;
				else {
					std::string	message = "Erreur dans cg_1to1_write step 3,  val2 = "+ std::to_string(val2) + " n'est ni au min ni au max " + std::to_string(isize2[ind2]);
					throw GMDSException (message);
				}

				transform[ind1] = ind2+1;
				if (is_val1_min != is_val2_min)
					transform[ind1] = -transform[ind1];

				cgsize_t pts[6] = {pts1[0],pts1[1],pts1[2],pts1[6],pts1[7],pts1[8]},
				         ptsdonnor[6] = {pts2[0],pts2[1],pts2[2],pts2[6],pts2[7],pts2[8]};


				if (cg_1to1_write(m_indexFile, m_indexBase, b.id() + 1, interface, zonename2, pts, ptsdonnor, transform, &index_tf) != CG_OK) {
					std::cout << cg_get_error() << std::endl;
				}
			}


			cgsize_t pts_bc[12];

			int type_bc = 0;
			int cpt_ff = 0;
			int cpt_par = 0;
			int cpt_axis = 0;
			int cpt_out = 0;
			for(auto n : face.getIDs<Node>()){
				if(farfield->value(n) == 1) cpt_ff++;
				if(paroi->value(n) == 1) cpt_par++;
				if (axis->value(n) == 1) cpt_axis++;
				if (out->value(n) == 1) cpt_out++;
			}
			if(cpt_ff == 4){
				type_bc = 1;
			}else if(cpt_par == 4){
				type_bc = 2;
			}else if(cpt_axis == 4){
				type_bc = 3;
			}else if(cpt_out == 4){
				type_bc = 4;
			}

			switch (iFace) {
			case 0:
				pts_bc[0] = 1;	// X1
				pts_bc[1] = 1;	// Y1
				pts_bc[2] = 1;	// Z1

				pts_bc[3] = discrI;	// X2
				pts_bc[4] = discrJ;	// Y2
				pts_bc[5] = 1;			// Z2
				break;
			case 2:
				pts_bc[0] = 1;
				pts_bc[1] = 1;
				pts_bc[2] = 1;

				pts_bc[3] = discrI;
				pts_bc[4] = 1;
				pts_bc[5] = discrK;
				break;
			case 5:
				pts_bc[0] = discrI;
				pts_bc[1] = 1;
				pts_bc[2] = 1;

				pts_bc[3] = discrI;
				pts_bc[4] = discrJ;
				pts_bc[5] = discrK;
				break;
			case 3:
				pts_bc[0] = 1;
				pts_bc[1] = discrJ;
				pts_bc[2] = 1;

				pts_bc[3] = discrI;
				pts_bc[4] = discrJ;
				pts_bc[5] = discrK;
				break;
			case 4:
				pts_bc[0] = 1;
				pts_bc[1] = 1;
				pts_bc[2] = 1;

				pts_bc[3] = 1;
				pts_bc[4] = discrJ;
				pts_bc[5] = discrK;
				break;
			case 1:
				pts_bc[0] = 1;	// X1
				pts_bc[1] = 1;	// Y1
				pts_bc[2] = discrK;	// Z1

				pts_bc[3] = discrI;	// X2
				pts_bc[4] = discrJ;	// Y2
				pts_bc[5] = discrK;	// Z2
				break;
			default: break;
			}

			std::string bcType_s;
			if (type_bc == 1) {
				bcType_s = "FARFIELD";
			}
			else if (type_bc == 2) {
				bcType_s = "PAROI";
			}
			else if (type_bc == 3) {
				bcType_s = "SYMETRIE";
			}
			else if (type_bc == 4) {
				bcType_s = "SORTIE";
			}
			else{
				bcType_s = "ORFN";
			}

			// Family name
			char bc_type[32];
			strcpy(bc_type, bcType_s.c_str());
			writeBoundaryCondition(id_bc, pts_bc, b.id() + 1, bc_type, face.id()+1);
		}

		delete[] x_coords;
		delete[] y_coords;
		delete[] z_coords;
	}
}
/*----------------------------------------------------------------------------*/
void CGNSWriter3D::writeBoundaryCondition(int &num_bc, cgsize_t* pts, int id_zone, char ABCtype[32], int AEdgeID) const{

	std::string name = ABCtype;
	//std::cout<<"\t -> bc of type "<<name<<std::endl;

	char bc_name[32];
	std::stringstream bctype_ss;
	bctype_ss <<ABCtype <<" "<<std::setw(5)<<std::setfill('0')<<AEdgeID;
	std::string bctype_s = bctype_ss.str();
	strcpy(bc_name,bctype_s.c_str());

	if(cg_boco_write(m_indexFile,m_indexBase,id_zone,bc_name,BCTypeUserDefined,PointRange,2,pts,&num_bc)
	    != CG_OK) {
		std::cout<<cg_get_error()<<std::endl;
	}
	if(cg_goto(m_indexFile,m_indexBase,"Zone_t", id_zone, "ZoneBC_t", 1, "BC_t", num_bc, "end")
	    != CG_OK) {
		std::cout<<cg_get_error()<<std::endl;
	}
	if(cg_famname_write(ABCtype)
	    != CG_OK) {
		std::cout<<cg_get_error()<<std::endl;
	}
}
/*----------------------------------------------------------------------------*/
void CGNSWriter3D::write(const std::string &AInFileName, const std::string &AOutFileName, const std::string &AWorkingDir){

	//std::cout<<"Start reading"<<std::endl;
/*
	gmds::IGMeshIOService ioService(m_blocks);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.setDataOptions(gmds::N|gmds::R);
	vtkReader.read(AInFileName);*/

	//std::cout<<"End reading"<<std::endl;

	//MeshDoctor doc(m_blocks);
	//doc.buildFacesAndR2F();
	//doc.buildF2R(m_blocks->getModel());

	m_block_grid = m_blocks->newVariable<std::vector<TCellID>,GMDS_REGION>("block_grid");

	VarDiscrI = m_blocks->getVariable<int,GMDS_REGION>("discrI");
	VarDiscrJ = m_blocks->getVariable<int,GMDS_REGION>("discrJ");
	VarDiscrK = m_blocks->getVariable<int,GMDS_REGION>("discrK");


	std::vector<Variable<int>*> varBlocks;

	for (auto b : m_blocks->regions()) {
		Variable<int>* blockID = m_blocks->getVariable<int,GMDS_NODE>("block"+std::to_string(b));
		//varBlocks.push_back(blockID);
		//m_block_grid->value(b).resize(VarDiscrI->value(b)+VarDiscrJ->value(b)+VarDiscrK->value(b));
		std::vector<TCellID> nodes;
		nodes.reserve(VarDiscrI->value(b)*VarDiscrJ->value(b)*VarDiscrK->value(b));
		nodes.resize(VarDiscrI->value(b)*VarDiscrJ->value(b)*VarDiscrK->value(b));
		for(auto n : m_blocks->nodes()){
			if(blockID->value(n) > 0){
				int i = blockID->value(n)-1;
				nodes[i] = n;
			}
		}
		for(auto n : nodes){
			//std::cout<<n<<std::endl;
			m_block_grid->value(b).push_back(n);
		}
	}

	axis = m_blocks->getOrCreateVariable<int,GMDS_NODE>("Symetrie");

	initialize(AOutFileName,AWorkingDir);
	writeZones();
	finalize(AWorkingDir);
}
/*----------------------------------------------------------------------------*/
void CGNSWriter3D::finalize(const std::string &AWorkingDir) const{

	std::cout<<"========================================================="<<std::endl;
	std::cout<<"Finalize CGNS Writer"<<std::endl;
	std::cout<<"========================================================="<<std::endl;

	cg_close(m_indexFile);
}
/*----------------------------------------------------------------------------*/
void CGNSWriter3D::_getIndicesIdAndVal(const int* ipnts1, const int* ipnts2, bool* filtre,
                                  int &ind, int &val)
{
	for (int i=0; i<3; i++){
		if (!filtre[i] && (ipnts1[i] == ipnts2[i])){
			ind = i;
			val = ipnts1[i];
			return;
		}
	}
}