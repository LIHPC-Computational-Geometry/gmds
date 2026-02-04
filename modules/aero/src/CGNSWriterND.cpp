/*----------------------------------------------------------------------------*/
//
// Created by calderans on 29/06/2023.
//
/*----------------------------------------------------------------------------*/
#include "gmds/aero/CGNSWriterND.h"
//#include "gmds/ig/MeshDoctor.h"
//#include "gmds/io/IGMeshIOService.h"
//#include "gmds/io/VTKReader.h"
#include <fstream>
#include <iomanip>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
using namespace aero;

// this handles CGNS built without CGNS_ENABLE_SCOPING
#if CG_BUILD_SCOPE == 0
#  define CG_Structured Structured
#  define CG_RealDouble RealDouble
#  define CG_BCTypeUserDefined BCTypeUserDefined
#  define CG_PointRange PointRange
#  define CG_Unstructured Unstructured
#  define CG_TRI_3 TRI_3
#endif
/*----------------------------------------------------------------------------*/
CGNSWriterND::CGNSWriterND(Blocking3D *ABlocking, int ADim)
  :m_blocks(ABlocking),m_mesh(nullptr), m_cellDim(ADim), m_physdim(ADim)
{}
/*----------------------------------------------------------------------------*/
CGNSWriterND::CGNSWriterND(Mesh *AMesh, int ADim)
  :m_mesh(nullptr),m_blocks(AMesh), m_cellDim(ADim), m_physdim(ADim)
{}
/*----------------------------------------------------------------------------*/
CGNSWriterND::CGNSWriterND()
  :m_mesh(nullptr)
{
	m_blocks = new Mesh(MeshModel(DIM3 | N | F | R | N2F | N2R | F2N | F2R | R2N | R2F));
}
/*----------------------------------------------------------------------------*/
CGNSWriterND::~CGNSWriterND()
{}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::initialize(const std::string &AOutFileName, const std::string &dir){

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

	char basename[AOutFileName.length()+1];

	strcpy(basename, AOutFileName.c_str());

	cg_base_write(m_indexFile,basename,m_cellDim,m_physdim,&m_indexBase);

}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::writeZones()
{

	std::cout << "=========================================================" << std::endl;
	std::cout << "Start writing zones" << std::endl;
	std::cout << "=========================================================" << std::endl;

	int index_tf = 0;
	int id_bc = 0;

	std::vector<Variable<int>*> bc_vars;
	std::vector<Variable<int>*> zone_vars;

	if(m_cellDim == 2){
		for(auto var : m_blocks->getAllVariables(GMDS_EDGE)){
			if(var->getName() != "mark"){
				bc_vars.push_back(m_blocks->getVariable<int,GMDS_EDGE>(var->getName()));
			}
		}
		for(auto var : m_blocks->getAllVariables(GMDS_FACE)){
			if(var->getName() != "discrI" && var->getName() != "discrJ" && var->getName() != "discrK"
		    	&& var->getName() != "block_grid" && var->getName() != "mark"){
				zone_vars.push_back(m_blocks->getVariable<int,GMDS_FACE>(var->getName()));
			}
		}
	}
	else if(m_cellDim == 3){
		for(auto var : m_blocks->getAllVariables(GMDS_FACE)){
			if(var->getName() != "mark"){
				bc_vars.push_back(m_blocks->getVariable<int,GMDS_FACE>(var->getName()));
			}
		}
		for(auto var : m_blocks->getAllVariables(GMDS_REGION)){
			if(var->getName() != "discrI" && var->getName() != "discrJ" && var->getName() != "discrK"
			    && var->getName() != "block_grid" && var->getName() != "mark"){
				zone_vars.push_back(m_blocks->getVariable<int,GMDS_REGION>(var->getName()));
			}
		}
	}

	std::vector<int> blocks;
	if(m_cellDim == 2){
		for(auto fid : m_blocks->faces()) blocks.push_back(fid);
	}else if(m_cellDim == 3){
		for(auto rid : m_blocks->regions()) blocks.push_back(rid);
	}

	std::cout<<"Nb blocks "<<blocks.size()<<std::endl;

	for (auto bid : blocks) {

		std::cout << "Writing zone " << bid+1 << std::endl;

		// Le nom de la zone correspond au nom du bloc, dans notre cas on l'appel juste
		//  "Block_<Block_ID>"
		std::string zone_type;
		for(auto var : zone_vars){
			if(var->value(bid) == 1)
				zone_type = var->getName();
		}
		std::stringstream ss;
		ss << zone_type << " " << std::setw(5) << std::setfill('0') << bid+1;     // Pour le moment on utilise les faces du mesh comme bloc
		std::string name = ss.str();
		char zonename[32];
		strcpy(zonename, name.c_str());

		int discrI = VarDiscrI->value(bid);
		int discrJ = VarDiscrJ->value(bid);
		int discrK = m_cellDim == 3 ? VarDiscrK->value(bid) : 1;

		// Ici on va remplir les infos sur le nombre de sommets, mailles etc du bloc
		cgsize_t *zone_size;

		if(m_cellDim == 2){
			zone_size = new cgsize_t[6];
			zone_size[0] = discrI;
			zone_size[1] = discrJ;
			zone_size[2] = discrI - 1;
			zone_size[3] = discrJ - 1;
			zone_size[4] = 0;
			zone_size[5] = 0;
		}
		else if(m_cellDim == 3) {
			zone_size = new cgsize_t[9];
			zone_size[0] = discrI;
			zone_size[1] = discrJ;
			zone_size[2] = discrK;
			zone_size[3] = discrI - 1;
			zone_size[4] = discrJ - 1;
			zone_size[5] = discrK - 1;
			zone_size[6] = 0;
			zone_size[7] = 0;
			zone_size[8] = 0;
		}
		else {
			throw GMDSException("Dimension "+std::to_string(m_cellDim)+"D not supported in CGNS writer.");
		}

		if(cg_zone_write(m_indexFile, m_indexBase, zonename, zone_size, CG_Structured, &m_indexZone) != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}

		int size = discrI * discrJ * discrK;
		auto x_coords = new double[size];
		auto y_coords = new double[size];
		auto z_coords = new double[size];

		std::vector<TCellID> b_grid = m_block_grid->value(bid);


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

		cg_coord_write(m_indexFile, m_indexBase, m_indexZone, CG_RealDouble, coord_nameX, x_coords, &index_coord);
		cg_coord_write(m_indexFile, m_indexBase, m_indexZone, CG_RealDouble, coord_nameY, y_coords, &index_coord);
		if(m_cellDim == 3)
			cg_coord_write(m_indexFile, m_indexBase, m_indexZone, CG_RealDouble, coord_nameZ, z_coords, &index_coord);

		//Writing Boundary Condition
		if(m_cellDim == 2){
			Face b = m_blocks->get<Face>(bid);
			for(int iEdge = 0; iEdge<4; iEdge++){

				writeConnections2D(b, iEdge, index_tf, zone_vars);

				writeBoundaryCondition2D(id_bc, b, iEdge, bc_vars);
			}

		}
		else if(m_cellDim == 3){
			Region b = m_blocks->get<Region>(bid);

			for(int iFace = 0; iFace<6; iFace++){

				writeConnections3D(b, iFace, index_tf, zone_vars);

				writeBoundaryCondition3D(id_bc, b, iFace, bc_vars);
			}
		}

		delete[] x_coords;
		delete[] y_coords;
		delete[] z_coords;

		delete[] zone_size;
	}
}
void
CGNSWriterND::writeConnections3D(const Region& Ablock, int iFace, int& index_tf, const std::vector<Variable<int>*>& zone_vars) const{
	Face face = Ablock.get<Face>()[iFace];

	if(face.get<Region>().size() == 2) {

		int id_voisin = face.getIDs<Region>()[0] == Ablock.id() ? face.getIDs<Region>()[1] : face.getIDs<Region>()[0];
		int discrI = VarDiscrI->value(Ablock.id());
		int discrJ = VarDiscrJ->value(Ablock.id());
		int discrK = VarDiscrK->value(Ablock.id());

		char interface[33];
		std::stringstream interface_ss;
		interface_ss << "F" << std::setw(5) << std::setfill('0') << face.id() + 1 << " (B" << std::setw(5) << std::setfill('0') << Ablock.id() + 1 << " & B"
		             << std::setw(5) << std::setfill('0') << id_voisin + 1 << ")";
		std::string interface_s = interface_ss.str();
		strcpy(interface, interface_s.c_str());
		std::cout << "\t -> connection " << interface_s << std::endl;

		std::string zone_type2;
		for (auto var : zone_vars) {
			if (var->value(id_voisin) == 1) zone_type2 = var->getName();
		}
		std::stringstream name2_ss;
		name2_ss << zone_type2 << std::setw(5) << std::setfill('0') << id_voisin + 1;
		std::string name2 = name2_ss.str();     // Pour le moment on utilise les faces du mesh comme bloc
		char zonename2[32];
		strcpy(zonename2, name2.c_str());

		cgsize_t pts1[9], pts2[9];
		int idmin1, idinter1, idmax1;

		int discrI_vois = VarDiscrI->value(id_voisin), discrJ_vois = VarDiscrJ->value(id_voisin), discrK_vois = VarDiscrK->value(id_voisin);

		switch (iFace) {
		case 0:
			pts1[0] = 1;     // X1
			pts1[1] = 1;     // Y1
			pts1[2] = 1;     // Z1

			pts1[3] = discrI;     // X2
			pts1[4] = 1;          // Y2
			pts1[5] = 1;          // Z2

			pts1[6] = discrI;     // X3
			pts1[7] = discrJ;     // Y3
			pts1[8] = 1;          // Z3

			idmin1 = 0;
			idinter1 = 1;
			idmax1 = 2;
			break;
		case 1:
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
		case 2:
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
		case 5:
			pts1[0] = 1;          // X1
			pts1[1] = 1;          // Y1
			pts1[2] = discrK;     // Z1

			pts1[3] = discrI;     // X2
			pts1[4] = 1;          // Y2
			pts1[5] = discrK;     // Z2

			pts1[6] = discrI;     // X3
			pts1[7] = discrJ;     // Y3
			pts1[8] = discrK;     // Z3

			idmin1 = 4;
			idinter1 = 5;
			idmax1 = 6;
			break;
		default: break;
		}

		int id_glob_n0 = Ablock.getIDs<Node>()[idmin1];
		int id_glob_n1 = Ablock.getIDs<Node>()[idmax1];
		int id_glob_inter1 = Ablock.getIDs<Node>()[idinter1];

		Region voisin = m_blocks->get<Region>(id_voisin);

		std::vector<int> indices2;
		indices2.resize(3);

		std::vector<TCellID> nid_vois = voisin.getIDs<Node>();
		for (int i = 0; i < 8; i++) {
			if (nid_vois[i] == id_glob_n0) {
				indices2[0] = i;
			}
			else if (nid_vois[i] == id_glob_inter1) {
				indices2[1] = i;
			}
			else if (nid_vois[i] == id_glob_n1) {
				indices2[2] = i;
			}
		}

		for (int i = 0; i < 3; i++) {
			switch (indices2[i]) {
			case 0:
				pts2[0 + (3 * i)] = 1;     // X1
				pts2[1 + (3 * i)] = 1;     // Y1
				pts2[2 + (3 * i)] = 1;     // Z1
				break;
			case 1:
				pts2[0 + (3 * i)] = discrI_vois;
				pts2[1 + (3 * i)] = 1;
				pts2[2 + (3 * i)] = 1;
				break;
			case 2:
				pts2[0 + (3 * i)] = discrI_vois;
				pts2[1 + (3 * i)] = discrJ_vois;
				pts2[2 + (3 * i)] = 1;
				break;
			case 3:
				pts2[0 + (3 * i)] = 1;
				pts2[1 + (3 * i)] = discrJ_vois;
				pts2[2 + (3 * i)] = 1;
				break;
			case 4:
				pts2[0 + (3 * i)] = 1;
				pts2[1 + (3 * i)] = 1;
				pts2[2 + (3 * i)] = discrK_vois;
				break;
			case 5:
				pts2[0 + (3 * i)] = discrI_vois;
				pts2[1 + (3 * i)] = 1;
				pts2[2 + (3 * i)] = discrK_vois;
				break;
			case 6:
				pts2[0 + (3 * i)] = discrI_vois;
				pts2[1 + (3 * i)] = discrJ_vois;
				pts2[2 + (3 * i)] = discrK_vois;
				break;
			case 7:
				pts2[0 + (3 * i)] = 1;
				pts2[1 + (3 * i)] = discrJ_vois;
				pts2[2 + (3 * i)] = discrK_vois;
				break;
			default: break;
			}
		}


		int transform[3];

		bool filtre[3] = {false, false, false};
		bool filtre_vois[3] = {false, false, false};

		int isize1[3] = {discrI, discrJ, discrK};
		int isize2[3] = {discrI_vois, discrJ_vois, discrK_vois};

		int ind1, ind2;
		int val1, val2;
		int ipnts1[3], ipnts2[3];
		for (int k = 0; k < 3; k++) {
			ipnts1[k] = pts1[k];
			ipnts2[k] = pts1[6 + k];
		}
		_getIndicesIdAndVal(ipnts1, ipnts2, filtre, ind1, val1);
		filtre[ind1] = true;

		for (int k = 0; k < 3; k++) {
			ipnts1[k] = pts2[k];
			ipnts2[k] = pts2[6 + k];
		}
		_getIndicesIdAndVal(ipnts1, ipnts2, filtre_vois, ind2, val2);
		filtre_vois[ind2] = true;

		bool is_val1_min;
		if (val1 == 1)
			is_val1_min = true;
		else if (isize1[ind1] == val1)
			is_val1_min = false;
		else {
			std::string message =
			   "Erreur dans cg_1to1_write step 1,  val1 = " + std::to_string(val1) + " n'est ni au min ni au max " + std::to_string(isize1[ind1]);
			throw GMDSException(message);
		}

		// idem avec val2
		bool is_val2_min;
		if (val2 == 1)
			is_val2_min = true;
		else if (isize2[ind2] == val2)
			is_val2_min = false;
		else {
			std::string message =
			   "Erreur dans cg_1to1_write step 1,  val2 = " + std::to_string(val2) + " n'est ni au min ni au max " + std::to_string(isize2[ind2]);
			throw GMDSException(message);
		}

		transform[ind1] = ind2 + 1;
		if (is_val1_min == is_val2_min) transform[ind1] = -transform[ind1];

		for (int k = 0; k < 3; k++) {
			ipnts1[k] = pts1[k];
			ipnts2[k] = pts1[3 + k];
		}
		_getIndicesIdAndVal(ipnts1, ipnts2, filtre, ind1, val1);
		filtre[ind1] = true;

		for (int k = 0; k < 3; k++) {
			ipnts1[k] = pts2[k];
			ipnts2[k] = pts2[3 + k];
		}
		_getIndicesIdAndVal(ipnts1, ipnts2, filtre_vois, ind2, val2);
		filtre_vois[ind2] = true;

		if (val1 == 1)
			is_val1_min = true;
		else if (isize1[ind1] == val1)
			is_val1_min = false;
		else {
			std::string message =
			   "Erreur dans cg_1to1_write step 2,  val1 = " + std::to_string(val1) + " n'est ni au min ni au max " + std::to_string(isize1[ind1]);
			throw GMDSException(message);
		}

		// idem avec val2
		if (val2 == 1)
			is_val2_min = true;
		else if (isize2[ind2] == val2)
			is_val2_min = false;
		else {
			std::string message =
			   "Erreur dans cg_1to1_write step 2,  val2 = " + std::to_string(val2) + " n'est ni au min ni au max " + std::to_string(isize2[ind2]);
			throw GMDSException(message);
		}

		transform[ind1] = ind2 + 1;
		if (is_val1_min != is_val2_min) transform[ind1] = -transform[ind1];

		for (int k = 0; k < 3; k++) {
			ipnts1[k] = pts1[3 + k];
			ipnts2[k] = pts1[6 + k];
		}
		_getIndicesIdAndVal(ipnts1, ipnts2, filtre, ind1, val1);
		filtre[ind1] = true;

		for (int k = 0; k < 3; k++) {
			ipnts1[k] = pts2[3 + k];
			ipnts2[k] = pts2[6 + k];
		}
		_getIndicesIdAndVal(ipnts1, ipnts2, filtre_vois, ind2, val2);
		filtre_vois[ind2] = true;

		if (val1 == 1)
			is_val1_min = true;
		else if (isize1[ind1] == val1)
			is_val1_min = false;
		else {
			std::string message =
			   "Erreur dans cg_1to1_write step 3,  val1 = " + std::to_string(val1) + " n'est ni au min ni au max " + std::to_string(isize1[ind1]);
			throw GMDSException(message);
		}

		// idem avec val2
		if (val2 == 1)
			is_val2_min = true;
		else if (isize2[ind2] == val2)
			is_val2_min = false;
		else {
			std::string message =
			   "Erreur dans cg_1to1_write step 3,  val2 = " + std::to_string(val2) + " n'est ni au min ni au max " + std::to_string(isize2[ind2]);
			throw GMDSException(message);
		}

		transform[ind1] = ind2 + 1;
		if (is_val1_min != is_val2_min) transform[ind1] = -transform[ind1];

		cgsize_t pts[6] = {pts1[0], pts1[1], pts1[2], pts1[6], pts1[7], pts1[8]}, ptsdonnor[6] = {pts2[0], pts2[1], pts2[2], pts2[6], pts2[7], pts2[8]};

		if (cg_1to1_write(m_indexFile, m_indexBase, Ablock.id() + 1, interface, zonename2, pts, ptsdonnor, transform, &index_tf) != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
	}
}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::writeConnections2D(const Face& Ablock, int iEdge, int& index_tf, const std::vector<Variable<int>*>& zone_vars) const
{

	Edge edge = Ablock.get<Edge>()[iEdge];

	if(edge.get<Face>().size() == 2) {

		int id_voisin = edge.getIDs<Face>()[0] == Ablock.id() ? edge.getIDs<Face>()[1] : edge.getIDs<Face>()[0];
		int discrI = VarDiscrI->value(Ablock.id());
		int discrJ = VarDiscrJ->value(Ablock.id());

		char interface[33];
		std::stringstream interface_ss;
		interface_ss << "F" << std::setw(5) << std::setfill('0') << edge.id() + 1 << " (B" << std::setw(5) << std::setfill('0') << Ablock.id() + 1 << " & B"
		             << std::setw(5) << std::setfill('0') << id_voisin + 1 << ")";
		std::string interface_s = interface_ss.str();
		strcpy(interface, interface_s.c_str());
		std::cout << "\t -> connection " << interface_s << std::endl;

		std::string zone_type2;
		for (auto var : zone_vars) {
			if (var->value(id_voisin) == 1) zone_type2 = var->getName();
		}
		std::stringstream name2_ss;
		name2_ss << zone_type2 << std::setw(5) << std::setfill('0') << id_voisin + 1;
		std::string name2 = name2_ss.str();     // Pour le moment on utilise les faces du mesh comme bloc
		char zonename2[32];
		strcpy(zonename2, name2.c_str());

		cgsize_t pts1[4], pts2[4];
		int idmin1, idmax1;

		int discrI_vois = VarDiscrI->value(id_voisin), discrJ_vois = VarDiscrJ->value(id_voisin);

		switch (iEdge) {
		case 0:
			pts1[0] = 1;
			pts1[1] = 1;
			pts1[2] = discrI;
			pts1[3] = 1;

			idmin1 = 0;
			idmax1 = 1;
			break;
		case 1:
			pts1[0] = discrI;
			pts1[1] = 1;
			pts1[2] = discrI;
			pts1[3] = discrJ;

			idmin1 = 1;
			idmax1 = 2;
			break;
		case 2:
			pts1[0] = 1;
			pts1[1] = discrJ;
			pts1[2] = discrI;
			pts1[3] = discrJ;

			idmin1 = 3;
			idmax1 = 2;
			break;
		case 3:
			pts1[0] = 1;
			pts1[1] = 1;
			pts1[2] = 1;
			pts1[3] = discrJ;

			idmin1 = 0;
			idmax1 = 3;
			break;
		}

		int id_glob_n0 = Ablock.getIDs<Node>()[idmin1];
		int id_glob_n1 = Ablock.getIDs<Node>()[idmax1];

		Face voisin = m_blocks->get<Face>(id_voisin);

		std::vector<int> indices2;
		indices2.resize(2);

		std::vector<TCellID> nid_vois = voisin.getIDs<Node>();
		for (int i = 0; i < 4; i++) {
			if (nid_vois[i] == id_glob_n0) {
				indices2[0] = i;
			}
			else if (nid_vois[i] == id_glob_n1) {
				indices2[1] = i;
			}
		}

		for (int i = 0; i < 2; i++) {
			switch (indices2[i]) {
			case 0:
				pts2[0 + (2 * i)] = 1;     // X1
				pts2[1 + (2 * i)] = 1;     // Y1
				break;
			case 1:
				pts2[0 + (2 * i)] = discrI_vois;
				pts2[1 + (2 * i)] = 1;
				break;
			case 2:
				pts2[0 + (2 * i)] = discrI_vois;
				pts2[1 + (2 * i)] = discrJ_vois;
				break;
			case 3:
				pts2[0 + (2 * i)] = 1;
				pts2[1 + (2 * i)] = discrJ_vois;
				break;
			default: break;
			}
		}

		int transform[2];

		int ind1, ind2;
		int val1, val2;
		for (int i = 0; i < 2; ++i) {
			if(pts1[i] == pts1[i+2]){
				val1 = pts1[i];
				for (int i2 = 0; i2 < 2; ++i2) {
					if(pts2[i2] == pts2[i2+2]){
						transform[i] = i2 + 1;
						val2 = pts2[i2];
						if(val1 == val2){
							transform[i] = -transform[i];
						}
					}
				}
			}
			else{
				val1 = pts1[i];
				for (int i2 = 0; i2 < 2; ++i2) {
					if(pts2[i2] != pts2[i2+2]){
						transform[i] = i2 + 1;
						val2 = pts2[i2];
						if(val1 != val2){
							transform[i] = -transform[i];
						}
					}
				}
			}
		}

		cgsize_t pts[4] = {pts1[0], pts1[1], pts1[2], pts1[3]}, ptsdonnor[4] = {pts2[0], pts2[1], pts2[2], pts2[3]};

		if (cg_1to1_write(m_indexFile, m_indexBase, Ablock.id() + 1, interface, zonename2, pts, ptsdonnor, transform, &index_tf) != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
	}
}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::writeBoundaryCondition3D(int &num_bc, const Region& Ablock, int iFace, const std::vector<Variable<int>*>& bc_vars) const{


	Face face = Ablock.get<Face>()[iFace];

	int discrI = VarDiscrI->value(Ablock.id());
	int discrJ = VarDiscrJ->value(Ablock.id());
	int discrK = VarDiscrK->value(Ablock.id());

	cgsize_t pts_bc[12];

	int type_bc = -1;
	for(int ivar = 0; ivar < bc_vars.size(); ivar++){
		if(bc_vars[ivar]->value(face.id()) == 1)
			type_bc = ivar;
	}
	if(type_bc != -1) {

		switch (iFace) {
		case 0:
			pts_bc[0] = 1;     // X1
			pts_bc[1] = 1;     // Y1
			pts_bc[2] = 1;     // Z1

			pts_bc[3] = discrI;     // X2
			pts_bc[4] = discrJ;     // Y2
			pts_bc[5] = 1;          // Z2
			break;
		case 1:
			pts_bc[0] = 1;
			pts_bc[1] = 1;
			pts_bc[2] = 1;

			pts_bc[3] = discrI;
			pts_bc[4] = 1;
			pts_bc[5] = discrK;
			break;
		case 2:
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
		case 5:
			pts_bc[0] = 1;          // X1
			pts_bc[1] = 1;          // Y1
			pts_bc[2] = discrK;     // Z1

			pts_bc[3] = discrI;     // X2
			pts_bc[4] = discrJ;     // Y2
			pts_bc[5] = discrK;     // Z2
			break;
		default: break;
		}

		char bc_type[32];
		std::string bcType_s;
		bcType_s = bc_vars[type_bc]->getName();
		strcpy(bc_type, bcType_s.c_str());

		char bc_name[32];
		std::stringstream bcname_ss;
		bcname_ss << bcType_s << " " << std::setw(5) << std::setfill('0') << face.id();
		std::string bcname_s = bcname_ss.str();
		strcpy(bc_name, bcname_s.c_str());

		if (cg_boco_write(m_indexFile, m_indexBase, Ablock.id() + 1, bc_name, CG_BCTypeUserDefined, CG_PointRange, 2, pts_bc, &num_bc) != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
		if (cg_goto(m_indexFile, m_indexBase, "Zone_t", Ablock.id() + 1, "ZoneBC_t", 1, "BC_t", num_bc, "end") != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
		if (cg_famname_write(bc_type) != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
	}
}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::writeBoundaryCondition2D(int &num_bc, const Face& Ablock, int iEdge, const std::vector<Variable<int>*>& bc_vars) const{
	Edge edge = Ablock.get<Edge>()[iEdge];

	int discrI = VarDiscrI->value(Ablock.id());
	int discrJ = VarDiscrJ->value(Ablock.id());

	cgsize_t pts_bc[4];

	int type_bc = -1;
	for(int ivar = 0; ivar < bc_vars.size(); ivar++){
		if(bc_vars[ivar]->value(edge.id()) == 1)
			type_bc = ivar;
	}

	if(type_bc != -1) {

		switch (iEdge) {
		case 0:
			pts_bc[0] = 1;
			pts_bc[1] = 1;
			pts_bc[2] = discrI;
			pts_bc[3] = 1;
			break;
		case 1:
			pts_bc[0] = discrI;
			pts_bc[1] = 1;
			pts_bc[2] = discrI;
			pts_bc[3] = discrJ;
			break;
		case 2:
			pts_bc[0] = 1;
			pts_bc[1] = discrJ;
			pts_bc[2] = discrI;
			pts_bc[3] = discrJ;
			break;
		case 3:
			pts_bc[0] = 1;
			pts_bc[1] = 1;
			pts_bc[2] = 1;
			pts_bc[3] = discrJ;
			break;
		default: break;
		}

		char bc_type[32];
		std::string bcType_s;
		bcType_s = bc_vars[type_bc]->getName();
		strcpy(bc_type, bcType_s.c_str());

		char bc_name[32];
		std::stringstream bcname_ss;
		bcname_ss << bcType_s << " " << std::setw(5) << std::setfill('0') << edge.id();
		std::string bcname_s = bcname_ss.str();
		strcpy(bc_name, bcname_s.c_str());

		if (cg_boco_write(m_indexFile, m_indexBase, Ablock.id() + 1, bc_name, CG_BCTypeUserDefined, CG_PointRange, 2, pts_bc, &num_bc) != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
		if (cg_goto(m_indexFile, m_indexBase, "Zone_t", Ablock.id() + 1, "ZoneBC_t", 1, "BC_t", num_bc, "end") != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
		if (cg_famname_write(bc_type) != CG_OK) {
			std::cout << cg_get_error() << std::endl;
		}
	}
}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::write(const std::string &AInFileName, const std::string &AOutFileName, const std::string &AWorkingDir){
	/*
	std::cout<<"Start reading"<<std::endl;

	gmds::IGMeshIOService ioService(m_blocks);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.setDataOptions(gmds::N|gmds::R);
	vtkReader.read(AInFileName);

	std::cout<<"End reading"<<std::endl;
*/
	//MeshDoctor doc(m_blocks);
	//doc.buildFacesAndR2F();
	//doc.buildF2R(m_blocks->getModel());

	if(m_cellDim == 2){
		m_block_grid = m_blocks->newVariable<std::vector<TCellID>,GMDS_FACE>("block_grid");
	}else if(m_cellDim == 3){
		m_block_grid = m_blocks->newVariable<std::vector<TCellID>,GMDS_REGION>("block_grid");
	}

	if(m_cellDim == 2){
		VarDiscrI = m_blocks->getVariable<int,GMDS_FACE>("discrI");
		VarDiscrJ = m_blocks->getVariable<int,GMDS_FACE>("discrJ");
	}else if(m_cellDim == 3){
		VarDiscrI = m_blocks->getVariable<int,GMDS_REGION>("discrI");
		VarDiscrJ = m_blocks->getVariable<int,GMDS_REGION>("discrJ");
		VarDiscrK = m_blocks->getVariable<int,GMDS_REGION>("discrK");
	}


	std::vector<Variable<int>*> varBlocks;


	std::vector<int> blocks;
	if(m_cellDim == 2){
		for(auto fid : m_blocks->faces()) blocks.push_back(fid);
	}else if(m_cellDim == 3){
		for(auto rid : m_blocks->regions()) blocks.push_back(rid);
	}

	for (auto b : blocks) {
		Variable<int>* blockID = m_blocks->getVariable<int,GMDS_NODE>("block"+std::to_string(b));
		std::vector<TCellID> nodes;
		int discrK = m_cellDim == 3 ? VarDiscrK->value(b) : 1;
		nodes.reserve(VarDiscrI->value(b)*VarDiscrJ->value(b)*discrK);
		nodes.resize(VarDiscrI->value(b)*VarDiscrJ->value(b)*discrK);
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

	initialize(AOutFileName,AWorkingDir);
	writeZones();
	finalize(AWorkingDir);
}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::finalize(const std::string &AWorkingDir) const{

	std::cout<<"========================================================="<<std::endl;
	std::cout<<"Finalize CGNS Writer"<<std::endl;
	std::cout<<"========================================================="<<std::endl;

	cg_close(m_indexFile);
}
/*----------------------------------------------------------------------------*/
void
CGNSWriterND::_getIndicesIdAndVal(const int* ipnts1, const int* ipnts2, bool* filtre,
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
