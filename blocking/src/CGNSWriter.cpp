/*----------------------------------------------------------------------------*/
//
// Created by calderans on 22/03/2022.
//
/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CGNSWriter.h"
#include <fstream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
using namespace blocking;

#define Structured CG_Structured
#define RealDouble CG_RealDouble
#define BCTypeUserDefined CG_BCTypeUserDefined
#define PointRange CG_PointRange
#define Unstructured CG_Unstructured
#define TRI_3 CG_TRI_3
/*----------------------------------------------------------------------------*/
CGNSWriter::CGNSWriter(Blocking2D *ABlocking)
			:m_blocks(ABlocking),m_mesh(nullptr)
{}
/*----------------------------------------------------------------------------*/
CGNSWriter::CGNSWriter(Mesh *AMesh)
  :m_mesh(AMesh),m_blocks(nullptr)
{}
/*----------------------------------------------------------------------------*/
CGNSWriter::~CGNSWriter()
{}
/*----------------------------------------------------------------------------*/
void CGNSWriter::initialize(const std::string &AFileName){

	std::cout<<"========================================================="<<std::endl;
	std::cout<<"Initialize CGNS Writer"<<std::endl;
	std::cout<<"========================================================="<<std::endl;

	m_stream = new std::ofstream(AFileName, std::ios::out);
	if (!m_stream){
		std::string s ="Impossible to create a CGNS File: "+AFileName;
		throw GMDSException(s);
	}

	char path[21];
	strcpy(path,"/home/calderans/dev/");

	cg_set_path(path);

	int n = AFileName.length();

	// declaring character array
	char char_array[n + 1];

	// copying the contents of the
	// string to char array
	strcpy(char_array, AFileName.c_str());

	cg_open("aeropipeline.cgns",CG_MODE_WRITE,&m_indexFile);

	char basename[33];
	strcpy(basename, "aeropipeline.cgns");

	m_cellDim = 2;
	m_physdim = 3;

	cg_base_write(m_indexFile,basename,m_cellDim,m_physdim,&m_indexBase);

}
/*----------------------------------------------------------------------------*/
void CGNSWriter::writeZones(){

	std::cout<<"========================================================="<<std::endl;
	std::cout<<"Start writing zones"<<std::endl;
	std::cout<<"========================================================="<<std::endl;

	int index_tf = 0;

	int id_bc = 0;

	for(auto b : m_blocks->allBlocks()){

		std::cout<<"Writing zone "<<b.id()<<std::endl;

		//Le nom de la zone correspond au nom du bloc, dans notre cas on l'appel juste
		// "Block_<Block_ID>"
		std::string name = "Block_"+std::to_string(b.id());//Pour le moment on utilise les faces du mesh comme bloc
		char zonename[32];
		strcpy(zonename,name.c_str());

		int discrI = b.getNbDiscretizationI();
		int discrJ = b.getNbDiscretizationJ();

		// Ici on va remplir les infos sur le nombre de sommets, mailles etc du bloc
		cgsize_t zone_size[6] = {discrI,
		                         discrJ,
		                         discrI - 1,
		                         discrJ - 1,
		                         0, 0};

		cg_zone_write(m_indexFile,m_indexBase,zonename,zone_size,Structured,&m_indexZone);


		double x_coords[discrI*discrJ];
		double y_coords[discrI*discrJ];
		double z_coords[discrI*discrJ];

		int i_coord = 0;
		for(int j = 0; j<discrJ; j++){
			for(int i = 0; i<discrI; i++){
				x_coords[i_coord] = b(i,j).X();
				y_coords[i_coord] = b(i,j).Y();
				z_coords[i_coord] = 0;
				i_coord++;
			}
		}

		int index_coord;

		char coord_nameX[32];
		strcpy(coord_nameX,"CoordinateX");
		char coord_nameY[32];
		strcpy(coord_nameY,"CoordinateY");
		char coord_nameZ[32];
		strcpy(coord_nameZ,"CoordinateZ");


		cg_coord_write(m_indexFile,m_indexBase,m_indexZone,RealDouble, coord_nameX, x_coords, &index_coord);
		cg_coord_write(m_indexFile,m_indexBase,m_indexZone,RealDouble, coord_nameY, y_coords, &index_coord);
		cg_coord_write(m_indexFile,m_indexBase,m_indexZone,RealDouble, coord_nameZ, z_coords, &index_coord);

		b.computeT();

		for(auto const &connection : b.getT()){

			int id = connection.first;
			std::vector<int> T = connection.second;

			Blocking2D::Block voisin = m_blocks->block(id);

			char interface[33];
			std::string interface_s = std::to_string(b.id())+"_to_"+ std::to_string(id);
			strcpy(interface, interface_s.c_str());
			std::cout<<"\t -> connection "<<interface_s<<std::endl;

			std::string name2 = "Block_"+std::to_string(id);//Pour le moment on utilise les faces du mesh comme bloc
			char zonename2[32];
			strcpy(zonename2,name2.c_str());

			cgsize_t pts1[4],pts2[4];

			std::cout<<b.getInterfaceInfo()[id][0]<<std::endl;
			std::cout<<b.getInterfaceInfo()[id][1]<<std::endl;

			switch (b.getInterfaceInfo()[id][0]) {
				case 0:
					pts1[0] = 1;
					pts1[1] = 1;
					pts1[2] = discrI;
					pts1[3] = 1;
					break;
				case 1:
					pts1[0] = discrI;
					pts1[1] = 1;
					pts1[2] = discrI;
					pts1[3] = discrJ;
					break;
				case 2:
					pts1[0] = 1;
					pts1[1] = discrJ;
					pts1[2] = discrI;
					pts1[3] = discrJ;
					break;
				case 3:
					pts1[0] = 1;
					pts1[1] = 1;
					pts1[2] = 1;
					pts1[3] = discrJ;
					break;
			}

			switch (b.getInterfaceInfo()[id][1]) {
				case 0:
					pts2[0] = 1;
					pts2[1] = 1;
					pts2[2] = voisin.getNbDiscretizationI();
					pts2[3] = 1;
					break;
				case 1:
					pts2[0] = voisin.getNbDiscretizationI();
					pts2[1] = 1;
					pts2[2] = voisin.getNbDiscretizationI();
					pts2[3] = voisin.getNbDiscretizationJ();
					break;
				case 2:
					pts2[0] = 1;
					pts2[1] = voisin.getNbDiscretizationJ();
					pts2[2] = voisin.getNbDiscretizationI();
					pts2[3] = voisin.getNbDiscretizationJ();
					break;
				case 3:
					pts2[0] = 1;
					pts2[1] = 1;
					pts2[2] = 1;
					pts2[3] = voisin.getNbDiscretizationJ();
					break;
			}

			int transform[2];
			transform[0] = b.getT()[id][0];
			transform[1] = b.getT()[id][1];

			if(cg_1to1_write(m_indexFile,m_indexBase,b.id()+1,interface,zonename2,pts1,pts2,transform,&index_tf)
			    != CG_OK) {
				std::cout<<cg_get_error()<<std::endl;
			}

		}

		//get coords of the 2 block corners that define boundary condition (bc), all mesh nodes between will be in the bc

		bool farfield_var,vehicule_var;

		bool farfield_found = false;
		bool vehicule_found = false;
		cgsize_t pts_bc[4];

		//For now, we assume that a bc is on only one edge of a block
		Variable<int>* farfield;
		Variable<int>* vehicule;

		try {
			farfield = m_blocks->getVariable<int, GMDS_NODE>("Farfield");
			farfield_var = true;
		} catch(GMDSException e){
			farfield_var = false;
		}
		try {
			vehicule = m_blocks->getVariable<int, GMDS_NODE>("Vehicule");
			vehicule_var = true;
		}catch (GMDSException e){
			vehicule_var = false;
		}
		std::vector<int> indices;

		if(farfield_var) {
			for (int i = 0; i < 4; i++) {
				Node n = b.getNode(i);
				if (farfield->value(n.id()) == 1) {
					farfield_found = true;
					indices.push_back(i);
				}
			}
			if(indices.size() == 1) indices.push_back(indices[0]); //Only one node of the block is in the bc, this may be a conception error: NEED TO VERIFY
			if (farfield_found) {
				for (int i = 0; i < indices.size(); i++) {
					switch (indices[i]) {
					case 0:
						pts_bc[0 + (i * 2)] = 1;
						pts_bc[1 + (i * 2)] = 1;
						break;
					case 1:
						pts_bc[0 + (i * 2)] = discrI;
						pts_bc[1 + (i * 2)] = 1;
						break;
					case 2:
						pts_bc[0 + (i * 2)] = discrI;
						pts_bc[1 + (i * 2)] = discrJ;
						break;
					case 3:
						pts_bc[0 + (i * 2)] = 1;
						pts_bc[1 + (i * 2)] = discrJ;
						break;
					}
				}

				char bc_name[32];
				strcpy(bc_name, "Farfield");

				writeBoundaryCondition(id_bc, pts_bc, b.id() + 1, bc_name);
			}
		}

		indices.clear();

		if(vehicule_var) {
			for (int i = 0; i < 4; i++) {
				Node n = b.getNode(i);
				if (vehicule->value(n.id()) == 1) {
					vehicule_found = true;
					indices.push_back(i);
				}
			}
			if(indices.size() == 1) indices.push_back(indices[0]);//Only one node of the block is in the bc, this may be a conception error: NEED TO VERIFY
			if (vehicule_found) {
				for (int i = 0; i < indices.size(); i++) {
					switch (indices[i]) {
					case 0:
						pts_bc[0 + (i * 2)] = 1;
						pts_bc[1 + (i * 2)] = 1;
						break;
					case 1:
						pts_bc[0 + (i * 2)] = discrI;
						pts_bc[1 + (i * 2)] = 1;
						break;
					case 2:
						pts_bc[0 + (i * 2)] = discrI;
						pts_bc[1 + (i * 2)] = discrJ;
						break;
					case 3:
						pts_bc[0 + (i * 2)] = 1;
						pts_bc[1 + (i * 2)] = discrJ;
						break;
					}
				}

				char bc_name[32];
				strcpy(bc_name, "Vehicule");

				writeBoundaryCondition(id_bc, pts_bc, b.id() + 1, bc_name);
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void CGNSWriter::writeBoundaryCondition(int &num_bc, cgsize_t* pts, int id_zone, char ABCtype[32]){

	std::string name = ABCtype;
	std::cout<<"\t -> bc of type "<<name<<std::endl;

	char bc_name[32];
	strcpy(bc_name,("BC_on_ENT" + std::to_string(num_bc)).c_str());

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
void CGNSWriter::writeTri(){

	std::string name = "Zone_tri";
	char zonename[32];
	strcpy(zonename,name.c_str());

	int nbNodes = m_mesh->getNbNodes();
	int nbCells = m_mesh->getNbFaces();

	cgsize_t zone_size[3];
	zone_size[0] = nbNodes;
	zone_size[1] = nbCells;
	zone_size[2] = 0;

	if(cg_zone_write(m_indexFile,m_indexBase,zonename,zone_size,Unstructured,&m_indexZone) != CG_OK){
		std::cout<<"Zone error : "<<cg_get_error()<<std::endl;
	}

	double x_coords[nbNodes];
	double y_coords[nbNodes];
	double z_coords[nbNodes];

	for(int i = 0; i<nbNodes; i++){
			x_coords[i] = m_mesh->get<Node>(i).X();
			y_coords[i] = m_mesh->get<Node>(i).Y();
			z_coords[i] = 0;
	}


	int index_coord;

	char coord_nameX[32];
	strcpy(coord_nameX,"CoordinateX");
	char coord_nameY[32];
	strcpy(coord_nameY,"CoordinateY");
	char coord_nameZ[32];
	strcpy(coord_nameZ,"CoordinateZ");


	cg_coord_write(m_indexFile,m_indexBase,m_indexZone,RealDouble, coord_nameX, x_coords, &index_coord);
	cg_coord_write(m_indexFile,m_indexBase,m_indexZone,RealDouble, coord_nameY, y_coords, &index_coord);
	cg_coord_write(m_indexFile,m_indexBase,m_indexZone,RealDouble, coord_nameZ, z_coords, &index_coord);

	int indexSection = 0;

	char sectionName[32];
	strcpy(sectionName,"Triangles");

	int i = 0;

	cgsize_t elems[nbCells*3];
	for(auto t : m_mesh->faces()){
		Face tri = m_mesh->get<Face>(t);
		elems[i] = tri.getIDs<Node>()[0]+1;
		elems[i+1] = tri.getIDs<Node>()[1]+1;
		elems[i+2] = tri.getIDs<Node>()[2]+1;
		i+=3;
	}

	if(cg_section_write(m_indexFile,m_indexBase,m_indexZone,sectionName,TRI_3,0,nbNodes-1,0,elems,&indexSection) != CG_OK){
		std::cout<<"Section error : "<<cg_get_error()<<std::endl;
	}
}
/*----------------------------------------------------------------------------*/
void CGNSWriter::finalize(){

	std::cout<<"========================================================="<<std::endl;
	std::cout<<"Finalize CGNS Writer"<<std::endl;
	std::cout<<"========================================================="<<std::endl;

	char path[32];
	strcpy(path,"/home/dev/aeropipeline.cgns");

	cg_save_as(m_indexFile,path,CG_FILE_HDF5,0);

	cg_close(m_indexFile);
}
/*----------------------------------------------------------------------------*/
void CGNSWriter::write(const std::string &AFileName){
	initialize(AFileName);
	writeZones();
	finalize();
}
