/*----------------------------------------------------------------------------*/
#include "gmds/blocking/Blocking.h"
/*----------------------------------------------------------------------------*/
#include <fstream>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/utils/Exception.h>

#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
/*----------------------------------------------------------------------------*/
Blocking::Blocking() {}
/*----------------------------------------------------------------------------*/
Blocking::~Blocking() {}
/*----------------------------------------------------------------------------*/
Blocking::STATUS
Blocking::execute()
{
	return Blocking::SUCCESS;
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid2d()
{
	this->createGrid2d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 3,3);
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid3d()
{
	this->createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 3,3,3);
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid2d(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy)
{
	Dart_handle* dhs = new Dart_handle[ANx*ANy];

	for(int i=0; i<ANx; i++) {
		for(int j=0; j<ANy; j++) {

			double x0 = APmin.X() + ((double) i / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
			double y0 = APmin.Y() + ((double) j / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));
			double x1 = APmin.X() + ((double) (i+1) / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
			double y1 = APmin.Y() + ((double) (j+1) / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));

			Dart_handle dh =
				lcc_.make_quadrangle(LCC_3::Point(x0, y0, 0.),
											LCC_3::Point(x1, y0, 0.),
											LCC_3::Point(x1, y1, 0.),
											LCC_3::Point(x0, y1, 0.));

			dhs[i*ANy+j] = dh;
		}
	}


	for(int i=0; i<ANx; i++) {
		for (int j=0; j<ANy; j++) {
			Dart_handle dh = dhs[i*ANy+j];

			if(i != ANx-1) {
				Dart_handle dhi = dhs[(i+1)*ANy+j];
				lcc_.sew<2>(lcc_.alpha(dh, 0, 1), lcc_.alpha(dhi, 1));
			}
			if(j != ANy-1) {
				Dart_handle dhj = dhs[i*ANy+j+1];
				lcc_.sew<2>(lcc_.alpha(dh, 1, 0, 1), dhj);
			}
		}
	}
	//lcc_.display_darts(std::cout, true);
	delete[] dhs;

	if (!lcc_.is_valid()) {
		std::string s ="Blocking::createGrid lcc not valid";
		throw gmds::GMDSException(s);
	}
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid3d(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy, int ANz)
{
	Dart_handle* dhs = new Dart_handle[ANx*ANy*ANz];

	for(int i=0; i<ANx; i++) {
		for(int j=0; j<ANy; j++) {
			for(int k=0; k<ANz; k++) {

				double x0 = APmin.X() + ((double) i / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
				double y0 = APmin.Y() + ((double) j / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));
				double z0 = APmin.Z() + ((double) k / (double) ANz) * (APmax.Z() + (-1. * APmin.Z()));
				double x1 = APmin.X() + ((double) (i+1) / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
				double y1 = APmin.Y() + ((double) (j+1) / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));
				double z1 = APmin.Z() + ((double) (k+1) / (double) ANz) * (APmax.Z() + (-1. * APmin.Z()));

				Dart_handle dh =
				   lcc_.make_hexahedron(LCC_3::Point(x0, y0, z0),
												LCC_3::Point(x1, y0, z0),
				                        LCC_3::Point(x1, y1, z0),
												LCC_3::Point(x0, y1, z0),
												LCC_3::Point(x0, y1, z1),
												LCC_3::Point(x0, y0, z1),
												LCC_3::Point(x1, y0, z1),
												LCC_3::Point(x1, y1, z1));

				dhs[i*(ANy*ANz)+j*ANz+k] = dh;
			}
		}
	}

	for(int i=0; i<ANx; i++) {
		for (int j=0; j<ANy; j++) {
			for (int k=0; k<ANz; k++) {
				Dart_handle dh = dhs[i*(ANy*ANz)+j*ANz+k];

				if(i != ANx-1) {
					Dart_handle dhi = dhs[(i+1)*(ANy*ANz)+j*ANz+k];
					lcc_.sew<3>(lcc_.alpha(dh, 1, 0, 1, 2), lcc_.alpha(dhi, 2));
				}
				if(j != ANy-1) {
					Dart_handle dhj = dhs[i*(ANy*ANz)+(j+1)*ANz+k];
					lcc_.sew<3>(lcc_.alpha(dh, 2, 1, 0, 1, 2), dhj);
				}
				if(k != ANz-1) {
					Dart_handle dhk = dhs[i*(ANy*ANz)+j*ANz+(k+1)];
					lcc_.sew<3>(lcc_.alpha(dh, 0, 1, 2), lcc_.alpha(dhk, 1, 2));
				}
			}
		}
	}

	delete[] dhs;

	if (!lcc_.is_valid()) {
		std::string s ="Blocking::createGrid lcc not valid";
		throw gmds::GMDSException(s);
	}
}
/*----------------------------------------------------------------------------*/
void Blocking::createBlocks3dFromMesh(const Mesh &AMesh){

	//Map use to store faces created in the LCC
	std::map<std::tuple<int,int,int,int>, Dart_handle> FtoD;

	//Get all regions from the mesh, we suppose the mesh contains only hexahedra
	for(auto r : AMesh.regions()){
		Region region = AMesh.get<Region>(r);

		createNewHex(region.get<Node>(),FtoD);
	}

	if (!lcc_.is_valid()) {
		std::string s ="Blocking::createBlocksFromMesh lcc not valid";
		throw gmds::GMDSException(s);
	}

}
/*----------------------------------------------------------------------------*/
void Blocking::createNewHex(std::vector<Node> ANodes, std::map<std::tuple<int,int,int,int>, Dart_handle> &AFtoD){

	//We consider that in GMDS a hex is represented as
	/*      7----------5
	 *     /|         /|
	 *    / |        / |
	 *   4--|-------5  |
	 *   |  3-------|--2
	 *   | /		    | /
	 *   |/		    |/
	 *   0----------1
	 */

	Dart_handle d1 = createNewQuad(ANodes[0],ANodes[1],ANodes[2],ANodes[3],AFtoD);
	Dart_handle d2 = createNewQuad(ANodes[0],ANodes[1],ANodes[5],ANodes[4],AFtoD);
	Dart_handle d3 = createNewQuad(ANodes[1],ANodes[2],ANodes[6],ANodes[5],AFtoD);
	Dart_handle d4 = createNewQuad(ANodes[2],ANodes[3],ANodes[7],ANodes[6],AFtoD);
	Dart_handle d5 = createNewQuad(ANodes[3],ANodes[0],ANodes[4],ANodes[7],AFtoD);
	Dart_handle d6 = createNewQuad(ANodes[4],ANodes[5],ANodes[6],ANodes[7],AFtoD);

	//sewing all faces of the hex
	lcc_.sew<2>(d1,d2);
	lcc_.sew<2>(lcc_.alpha<0,1>(d1),d3);
	lcc_.sew<2>(lcc_.alpha<0,1,0,1>(d1),d4);
	lcc_.sew<2>(lcc_.alpha<1,0>(d1),d5);

	lcc_.sew<2>(lcc_.alpha<0,1>(d2),lcc_.alpha<1>(d3));
	lcc_.sew<2>(lcc_.alpha<1>(d2),lcc_.alpha<0,1>(d5));
	lcc_.sew<2>(lcc_.alpha<1,0,1>(d2),d6);

	lcc_.sew<2>(lcc_.alpha<0,1>(d3),lcc_.alpha<1>(d4));
	lcc_.sew<2>(lcc_.alpha<1,0,1>(d3),lcc_.alpha<0,1>(d6));

	lcc_.sew<2>(lcc_.alpha<0,1>(d4),lcc_.alpha<1>(d5));
	lcc_.sew<2>(lcc_.alpha<1,0,1>(d4),lcc_.alpha<0,1,0,1>(d6));

	lcc_.sew<2>(lcc_.alpha<1,0,1>(d5),lcc_.alpha<1,0>(d6));

}
/*----------------------------------------------------------------------------*/
Dart_handle Blocking::createNewQuad(const Node n0,const Node n1,const Node n2,const Node n3, std::map<std::tuple<int,int,int,int>, Dart_handle> &AFtoD){

	Dart_handle d = lcc_.make_quadrangle(LCC_3::Point(n0.point().X(), n0.point().Y(), n0.point().Z()),
	                                     LCC_3::Point(n1.point().X(), n1.point().Y(), n1.point().Z()),
	                                     LCC_3::Point(n2.point().X(), n2.point().Y(), n2.point().Z()),
	                                     LCC_3::Point(n3.point().X(), n3.point().Y(), n3.point().Z()));

	int n0_id = n0.id();
	int n1_id = n1.id();
	int n2_id = n2.id();
	int n3_id = n3.id();

	Dart_handle dmin;
	std::tuple<int,int,int,int> f, fopp;

	//We have created a face in the LCC that only represent one side of a GMDS::Face
	//We are going to order this face and store it in order to sew it with the other side face

	//We set a global origin of the face by getting the vertex of minimum id
	//and the opposite face which means the same face oriented in the opposite direction
	if(n0_id<n1_id){
		if(n0_id<n2_id){
			if(n0_id<n3_id){//n0 is the min
				if(n1_id<n3_id){
					f = {n0_id,n1_id,n2_id,n3_id};
					fopp = {n0_id,n3_id,n2_id,n1_id};
					dmin = d;
				}else{
					f = {n0_id,n3_id,n2_id,n1_id};
					fopp = {n0_id,n1_id,n2_id,n3_id};
					dmin = lcc_.alpha<1>(d);
				}
			}else{//n3 is the min
				if(n0_id<n2_id){
					f = {n3_id,n0_id,n1_id,n2_id};
					fopp = {n3_id,n2_id,n1_id,n0_id};
					dmin = lcc_.alpha<1,0>(d);
				}else{
					f = {n3_id,n2_id,n1_id,n0_id};
					fopp = {n3_id,n0_id,n1_id,n2_id};
					dmin = lcc_.alpha<1,0,1>(d);
				}
			}
		}else{
			if(n2_id<n3_id){//n2 is the min
				if(n1_id<n3_id){
					f = {n2_id,n1_id,n0_id,n3_id};
					fopp = {n2_id,n3_id,n0_id,n1_id};
					dmin = lcc_.alpha<0,1,0>(d);
				}else{
					f = {n2_id,n3_id,n0_id,n1_id};
					fopp = {n2_id,n1_id,n0_id,n3_id};
					dmin = lcc_.alpha<0,1,0,1>(d);
				}
			}else{//n3 is the min
				if(n0_id<n2_id){
					f = {n3_id,n0_id,n1_id,n2_id};
					fopp = {n3_id,n2_id,n1_id,n0_id};
					dmin = lcc_.alpha<1,0>(d);
				}else{
					f = {n3_id,n2_id,n1_id,n0_id};
					fopp = {n3_id,n0_id,n1_id,n2_id};
					dmin = lcc_.alpha<1,0,1>(d);
				}
			}
		}
	}else{
		if(n1_id<n2_id){
			if(n1_id<n3_id){//n1 is the min
				if(n0_id<n2_id){
					f = {n1_id,n0_id,n3_id,n2_id};
					fopp = {n1_id,n2_id,n3_id,n0_id};
					dmin = lcc_.alpha<0>(d);
				}else{
					f = {n1_id,n2_id,n3_id,n0_id};
					fopp = {n1_id,n0_id,n3_id,n2_id};
					dmin = lcc_.alpha<0,1>(d);
				}
			}else{//n3 is the min
				if(n0_id<n2_id){
					f = {n3_id,n0_id,n1_id,n2_id};
					fopp = {n3_id,n2_id,n1_id,n0_id};
					dmin = lcc_.alpha<1,0>(d);
				}else{
					f = {n3_id,n2_id,n1_id,n0_id};
					fopp = {n3_id,n0_id,n1_id,n2_id};
					dmin = lcc_.alpha<1,0,1>(d);
				}
			}
		}else{
			if(n2_id<n3_id){//n2 is the min
				if(n1_id<n3_id){
					f = {n2_id,n1_id,n0_id,n3_id};
					fopp = {n2_id,n3_id,n0_id,n1_id};
					dmin = lcc_.alpha<0,1,0>(d);
				}else{
					f = {n2_id,n3_id,n0_id,n1_id};
					fopp = {n2_id,n1_id,n0_id,n3_id};
					dmin = lcc_.alpha<0,1,0,1>(d);
				}
			}else{//n3 is the min
				if(n0_id<n2_id){
					f = {n3_id,n0_id,n1_id,n2_id};
					fopp = {n3_id,n2_id,n1_id,n0_id};
					dmin = lcc_.alpha<1,0>(d);
				}else{
					f = {n3_id,n2_id,n1_id,n0_id};
					fopp = {n3_id,n0_id,n1_id,n2_id};
					dmin = lcc_.alpha<1,0,1>(d);
				}
			}
		}
	}

	for(auto f_test : AFtoD){
		//Search if the other side of the face is already created
		if(f == f_test.first){
			//Same orientation
			lcc_.sew<3>(dmin, f_test.second);
		}else if(fopp == f_test.first){
			//Opposite direction
			lcc_.sew<3>(lcc_.alpha<1>(dmin), f_test.second);
		}
	}

	AFtoD[f] = dmin;

	return d;
}

/*----------------------------------------------------------------------------*/
void Blocking::readVTKFile(std::string AFileName)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::F | gmds::R | gmds::F2N | gmds::R2N));

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F | gmds::R);
	vtkReader.read(AFileName);

	// check validity
	// we only handle meshes containing either faces or regions, not both, at the moment
	// and the mesh must have at least one
	if (((m.getNbFaces() == 0) && (m.getNbRegions() == 0)) ||
	    ((m.getNbFaces() != 0) && (m.getNbRegions() != 0))) {
		std::string s ="Blocking::readVTKFile we only handle meshes containing either faces or regions, "
		                "not both, at the moment and the mesh must have at least one";
		throw gmds::GMDSException(s);
	}

	createBlocks3dFromMesh(m);
}
/*----------------------------------------------------------------------------*/
void Blocking::writeMokaFile(std::string AFileName) const
{
	std::ofstream stream(AFileName, std::ios::out);
	if (!stream.is_open()){
		std::string s ="Impossible to create a Moka File (ASCII format): "+AFileName;
		throw gmds::GMDSException(s);
	}

	stream << "Moka file [ascii]" << std::endl;;

	// TODO investigate the magic_number
	const int magic_number = 128;
	stream << magic_number <<" 7 0 0 0 0 0 0" << std::endl;

	for(auto it = lcc_.darts().begin(); it != lcc_.darts().end(); ++it) {

		stream << lcc_.darts().index(lcc_.get_alpha<0>(it)) << " ";
		stream << lcc_.darts().index(lcc_.get_alpha<1>(it)) << " ";
		stream << lcc_.darts().index(lcc_.get_alpha<2>(it)) << " ";
		stream << lcc_.darts().index(lcc_.get_alpha<3>(it)) << " ";

		stream << magic_number << " 0 0 0 ";

		// check whether the dart is associated to a point
		LCC_3::Vertex_attribute_const_handle v = lcc_.vertex_attribute(it);
		if (v->dart() == it) {
			stream << "1 " << lcc_.point_of_vertex_attribute(v) << std::endl;
		}
		else {
			stream << "0" << std::endl;
		}
	}
	stream.close();
}
/*----------------------------------------------------------------------------*/
void Blocking::writeVTKFile(std::string AFileName) const
{
	// validity checks
	if (!lcc_.is_valid()) {
		std::string s("Blocking::writeVTKFile lcc is not valid");
		throw gmds::GMDSException(s);
	}

	gmds::MeshModel model(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N));
	gmds::Mesh m(model);

	std::map<LCC_3::Vertex_attribute_range::const_iterator, gmds::TCellID> v2n;

	// Display all the vertices of the map.
	for (LCC_3::Vertex_attribute_range::const_iterator
			  it=lcc_.vertex_attributes().begin(),
			  itend=lcc_.vertex_attributes().end();
		  it!=itend; ++it) {
		gmds::math::Point pt(lcc_.point_of_vertex_attribute(it).x(), lcc_.point_of_vertex_attribute(it).y(), lcc_.point_of_vertex_attribute(it).z());
		v2n[it] = m.newNode(pt).id();
	}

	for (auto it  = lcc_.one_dart_per_cell<2>().begin(); it != lcc_.one_dart_per_cell<2>().end(); it++) {
//		std::cout << "poyop " << lcc_.darts_of_orbit<0, 1>(it).size() << std::endl;

		int nb = lcc_.darts_of_orbit<0, 1>(it).size();
		switch (nb) {
		case 6: { // triangle
			LCC_3::Dart_const_handle d0 = it;
			LCC_3::Dart_const_handle d1 = lcc_.alpha(d0, 0, 1);
			LCC_3::Dart_const_handle d2 = lcc_.alpha(d1, 0, 1);
			LCC_3::Vertex_attribute_const_handle v0 = lcc_.vertex_attribute(d0);
			LCC_3::Vertex_attribute_const_handle v1 = lcc_.vertex_attribute(d1);
			LCC_3::Vertex_attribute_const_handle v2 = lcc_.vertex_attribute(d2);
			m.newTriangle(v2n[v0], v2n[v1], v2n[v2]);
		} break;
		case 8: { //quadrangle
			LCC_3::Dart_const_handle d0 = it;
			LCC_3::Dart_const_handle d1 = lcc_.alpha(d0, 0, 1);
			LCC_3::Dart_const_handle d2 = lcc_.alpha(d1, 0, 1);
			LCC_3::Dart_const_handle d3 = lcc_.alpha(d2, 0, 1);
			LCC_3::Vertex_attribute_const_handle v0 = lcc_.vertex_attribute(d0);
			LCC_3::Vertex_attribute_const_handle v1 = lcc_.vertex_attribute(d1);
			LCC_3::Vertex_attribute_const_handle v2 = lcc_.vertex_attribute(d2);
			LCC_3::Vertex_attribute_const_handle v3 = lcc_.vertex_attribute(d3);
			m.newQuad(v2n[v0], v2n[v1], v2n[v2], v2n[v3]);
		} break;
		default: {
			std::string s = "Blocking::writeVTKFile 2d cell has an unexpected number of darts " + std::to_string(nb);
			throw gmds::GMDSException(s);
		}
		}
	}

	// In some (all?) cases when in 2d a unique 3d cell exists and all darts self reference themselves by alpha3
	// check whether we have "real" 3d cells that we need to effectively write
	bool in2d = true;
	if(lcc_.one_dart_per_cell<3>().size() == 1) {
		int index = 0;
		for(auto it = lcc_.darts().begin(); it != lcc_.darts().end(); ++it) {
			if(index != lcc_.darts().index(lcc_.get_alpha<3>(it))) {
				in2d = false;
				break;
			}
			index++;
		}
	} else {
		in2d = false;
	}

	if(!in2d) {

		for (auto it = lcc_.one_dart_per_cell<3>().begin(); it != lcc_.one_dart_per_cell<3>().end(); it++) {
			int nb = lcc_.darts_of_orbit<0, 1, 2>(it).size();
			switch (nb) {
			case 24: {     // tetrahedron
				LCC_3::Dart_const_handle d0 = it;
				LCC_3::Dart_const_handle d1 = lcc_.alpha(d0, 0, 1);
				LCC_3::Dart_const_handle d2 = lcc_.alpha(d1, 0, 1);
				LCC_3::Dart_const_handle d3 = lcc_.alpha(d0, 2, 1, 0);
				LCC_3::Vertex_attribute_const_handle v0 = lcc_.vertex_attribute(d0);
				LCC_3::Vertex_attribute_const_handle v1 = lcc_.vertex_attribute(d1);
				LCC_3::Vertex_attribute_const_handle v2 = lcc_.vertex_attribute(d2);
				LCC_3::Vertex_attribute_const_handle v3 = lcc_.vertex_attribute(d3);
				m.newTet(v2n[v0], v2n[v1], v2n[v2], v2n[v3]);
			} break;
			case 36: {     // prism3
				            //			LCC_3::Dart_const_handle d0 = it;
				            //			LCC_3::Dart_const_handle d1 = lcc_.alpha(d0, 0, 1);
				            //			LCC_3::Dart_const_handle d2 = lcc_.alpha(d1, 0, 1);
				            //			LCC_3::Dart_const_handle d3 = lcc_.alpha(d0, 2, 1, 0);
				            //			LCC_3::Vertex_attribute_const_handle v0 = lcc_.vertex_attribute(d0);
				            //			LCC_3::Vertex_attribute_const_handle v1 = lcc_.vertex_attribute(d1);
				            //			LCC_3::Vertex_attribute_const_handle v2 = lcc_.vertex_attribute(d2);
				            //			LCC_3::Vertex_attribute_const_handle v3 = lcc_.vertex_attribute(d3);
				            //			m.newTet(v2n[v0], v2n[v1], v2n[v2], v2n[v3]);
			} break;
			case 48: {     // hexahedron
				LCC_3::Dart_const_handle d0 = it;
				LCC_3::Dart_const_handle d1 = lcc_.alpha(d0, 0, 1);
				LCC_3::Dart_const_handle d2 = lcc_.alpha(d1, 0, 1);
				LCC_3::Dart_const_handle d3 = lcc_.alpha(d2, 0, 1);
				LCC_3::Dart_const_handle d4 = lcc_.alpha(it, 2, 1, 0, 1, 2);
				LCC_3::Dart_const_handle d5 = lcc_.alpha(d4, 0, 1);
				LCC_3::Dart_const_handle d6 = lcc_.alpha(d5, 0, 1);
				LCC_3::Dart_const_handle d7 = lcc_.alpha(d6, 0, 1);
				LCC_3::Vertex_attribute_const_handle v0 = lcc_.vertex_attribute(d0);
				LCC_3::Vertex_attribute_const_handle v1 = lcc_.vertex_attribute(d1);
				LCC_3::Vertex_attribute_const_handle v2 = lcc_.vertex_attribute(d2);
				LCC_3::Vertex_attribute_const_handle v3 = lcc_.vertex_attribute(d3);
				LCC_3::Vertex_attribute_const_handle v4 = lcc_.vertex_attribute(d4);
				LCC_3::Vertex_attribute_const_handle v5 = lcc_.vertex_attribute(d5);
				LCC_3::Vertex_attribute_const_handle v6 = lcc_.vertex_attribute(d6);
				LCC_3::Vertex_attribute_const_handle v7 = lcc_.vertex_attribute(d7);
				m.newHex(v2n[v0], v2n[v1], v2n[v2], v2n[v3], v2n[v4], v2n[v5], v2n[v6], v2n[v7]);
			} break;
			default: {
				std::string s = "Blocking::writeVTKFile 3d cell has an unexpected number of darts " + std::to_string(nb);
				throw gmds::GMDSException(s);
			}
			}
		}
	}

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	if(!in2d) {
		vtkWriter.setCellOptions(gmds::N | gmds::R);
	} else {
		vtkWriter.setCellOptions(gmds::N|gmds::F);
	}
	vtkWriter.write(AFileName);
}
/*----------------------------------------------------------------------------*/
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/