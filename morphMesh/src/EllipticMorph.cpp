//
// Created by calderans on 30/08/23.
//

#include "gmds/morphMesh/EllipticMorph.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/morphMesh/FastLocalize.h"
#include <fstream>
#include <regex>


#ifdef WITH_LIMA
	#include "gmds/io/LimaReader.h"
	#include "gmds/io/LimaWriter.h"
#	include "gmds/math/Ray.h"
#	include "gmds/math/Triangle.h"
#endif //WITH_LIMA


/*----------------------------------------------------------------------------*/

using namespace gmds;
using  namespace morphmesh;

/*----------------------------------------------------------------------------*/
EllipticMorph::EllipticMorph(const std::string& AFilename, Mesh* AMesh):
  m_mesh(AMesh){

	   auto m_stream= new std::ifstream(AFilename.c_str(),std::ios::in);
	   if(!m_stream->is_open()){
		   std::string mess ="Impossible to read file "+AFilename;
		   throw GMDSException(mess);
	   }

	   std::string line;
	   bool found = false;
	   while (!found && std::getline(*m_stream,line))  {
		   found = std::regex_match(line, std::regex("Input file name"));
	   }
	   std::getline(*m_stream,line);
	   m_inputName = line;

	   //-----------------------------------------------------------------------------------
	   found = false;
	   while (!found && std::getline(*m_stream,line))  {
		   found = std::regex_match(line, std::regex("Output file name"));
	   }
	   std::getline(*m_stream,line);
	   m_outputName = line;

	   //-----------------------------------------------------------------------------------
	   found = false;
	   while (!found && std::getline(*m_stream,line))  {
		   found = std::regex_match(line, std::regex("To morph"));
	   }
	   while (std::getline(*m_stream,line) && !std::regex_match(line, std::regex("")))  {
		   m_to_morph.push_back(line);
	   }

	   //-----------------------------------------------------------------------------------
	   found = false;
	   while (!found && std::getline(*m_stream,line))  {
		   found = std::regex_match(line, std::regex("External lock"));
	   }
	   while (std::getline(*m_stream,line) && !std::regex_match(line, std::regex("")))  {
		   m_ext_lock.push_back(line);
	   }

	   //-----------------------------------------------------------------------------------
	   found = false;
	   while (!found && std::getline(*m_stream,line))  {
		   found = std::regex_match(line, std::regex("Internal lock"));
	   }
	   while (std::getline(*m_stream,line) && !std::regex_match(line, std::regex("")))  {
		   m_int_lock.push_back(line);
	   }

	   //-----------------------------------------------------------------------------------
	   found = false;
	   while (!found && std::getline(*m_stream,line))  {
		   found = std::regex_match(line, std::regex("Ellipses"));
	   }
	   std::string current_word;
	   *m_stream >> current_word;
	   while (!m_stream->eof()){

		   double X, coefa, coefb, coefc;

		   X = std::stod(current_word);
		   *m_stream >> current_word;
		   coefa = std::stod(current_word);
		   *m_stream >> current_word;
		   coefb = std::stod(current_word);
		   *m_stream >> current_word;
		   coefc = std::stod(current_word);

		   m_ellipses.push_back({X,coefa,coefb,coefc});

		   *m_stream >> current_word;
	   }
	   if(m_inputName.find(".vtk") != -1) {
		   std::cout<<"\t -> Read VTK format"<<std::endl;
		   gmds::IGMeshIOService ioService(m_mesh);
		   gmds::VTKReader vtkReader(&ioService);
		   vtkReader.setCellOptions(gmds::N | gmds::R);
		   vtkReader.setDataOptions(gmds::N | gmds::R);
		   vtkReader.read(m_inputName);
	   }
	   else if(m_inputName.find(".mli") != -1) {
#ifdef WITH_LIMA
		   std::cout<<"\t -> Read Lima format"<<std::endl;
		   gmds::LimaReader limaReader(*m_mesh);
		   limaReader.read(m_inputName, R|N);
#endif // WITH_LIMA
	   }else{
		   throw GMDSException("File format not supported, only .vtk or .mli2 supported");
	   }

	   gmds::MeshDoctor doc(m_mesh);
	   doc.buildN2R(m_mesh->getModel());
	   doc.buildFacesAndR2F();
	   doc.buildN2F(m_mesh->getModel());
	   doc.buildE();
	   doc.buildN2E(m_mesh->getModel());
	   doc.updateUpwardConnectivity();

	   m_locked = m_mesh->newMark<Node>();
	   m_morphRegions = m_mesh->newMark<Region>();
	   m_lockRegions = m_mesh->newMark<Region>();
}
/*----------------------------------------------------------------------------*/
EllipticMorph::~EllipticMorph() = default;
/*----------------------------------------------------------------------------*/
void EllipticMorph::execute()
{

	std::cout<<"============================== Starting Elliptic Morphing =============================="<<std::endl;

	markLockedCells();

	std::vector<TCellID> modified_nodes;

	if(m_ext_lock.empty()){
		modified_nodes = noExteriorLock();
	}else{
		modified_nodes = withExteriorLock();
	}

	for(auto n : modified_nodes){

		math::Vector3d vec(m_vecs[n]);
		Node node = m_mesh->get<Node>(n);
		double y = node.Y();
		double z = node.Z();

		node.setY(y+vec.Y());
		node.setZ(z+vec.Z());
	}

	finalize();

	std::cout<<"============================== Elliptic Morphing Finished =============================="<<std::endl;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> EllipticMorph::withExteriorLock(){

	FastLocalize fl_int(m_lockedIn);
	FastLocalize fl_ext(m_lockedOut);
	FastLocalize fl_morphed(m_morphed);

	for(auto f : m_mesh->faces()){
		Face face = m_mesh->get<Face>(f);
		if(face.getIDs<Region>().size() == 2){
			if ((m_mesh->isMarked<Region>(face.getIDs<Region>()[0],m_lockRegions) && !m_mesh->isMarked<Region>(face.getIDs<Region>()[1],m_lockRegions))
			    || (!m_mesh->isMarked<Region>(face.getIDs<Region>()[0],m_lockRegions) && m_mesh->isMarked<Region>(face.getIDs<Region>()[1],m_lockRegions))){
				outerSkin.push_back(face);
			}
		}
	}

	std::set<TCellID> internal_modified_nodes;

	for(int i = 0; i < m_ellipses.size()-1; i++) {

		double Xmin = m_ellipses[i][0];
		double Xmax = m_ellipses[i + 1][0];

		std::vector<Node> interval_morphed;

		for (auto n : m_morphed) {
			if (n.X() >= Xmin && n.X() < Xmax) {
				interval_morphed.push_back(n);
			}
		}

		double coefa_min = m_ellipses[i][1];
		double coefb_min = m_ellipses[i][2];
		double coefc_min = m_ellipses[i][3];

		double coefa_max = m_ellipses[i + 1][1];
		double coefb_max = m_ellipses[i + 1][2];
		double coefc_max = m_ellipses[i + 1][3];

		for (auto n : interval_morphed) {
			math::Point p = n.point();

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();

			double distXmin = p.X() - Xmin;
			double distXmax = Xmax - Xmin;
			double w = distXmin / distXmax;

			double coefy_min = p.Y() >= 0 ? coefa_min : coefc_min;
			double coefy_max = p.Y() >= 0 ? coefa_max : coefc_max;

			double coefy_test = coefy_min * (3 * pow(1 - w, 2) - 2 * pow(1 - w, 3)) + coefy_max * (3 * pow(w, 2) - 2 * pow(w, 3));
			double coefz_test = coefb_min * (3 * pow(1 - w, 2) - 2 * pow(1 - w, 3)) + coefb_max * (3 * pow(w, 2) - 2 * pow(w, 3));

			double theta = atan(n.Y() / n.Z());
			if (theta < 0) {
				theta *= -1;
			}

			double sinus = sin(theta);
			double cosinus = cos(theta);

			double newY = (p.Y() * ((coefy_test * sinus) + coefz_test * (1 - sinus)));
			double newZ = (p.Z() * ((coefz_test * cosinus) + coefy_test * (1 - cosinus)));

			math::Vector3d vec;
			vec.setXYZ(0, newY - p.Y(), newZ - p.Z());
			m_vecs[n.id()] = vec;
			internal_modified_nodes.insert(n.id());
		}
	}

	for(int i = 0; i < m_ellipses.size()-1; i++) {

		double Xmin = m_ellipses[i][0];
		double Xmax = m_ellipses[i+1][0];

		std::vector<Node> execution;

		for(auto n : m_mesh->nodes()){

			Node node = m_mesh->get<Node>(n);

			if (node.X() >= Xmin && node.X() < Xmax && !m_mesh->isMarked<Node>(n, m_locked)) {
				execution.push_back(node);
			}
		}

		for (const auto& n : execution) {
			internal_modified_nodes.insert(n.id());

			math::Point p = n.point();

			Node nearest_morphed = m_mesh->get<Node>(fl_morphed.find(p));

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			math::Point p_ext = m_mesh->get<Node>(fl_ext.find(p)).point();
			math::Point p_m = nearest_morphed.point();
			math::Point p_axe(p.X(), 0, 0);

			math::Ray ray(p_axe, p);

			Face nearest_locked_face;
			Edge nearest_locked_edge;
			math::Point p_intersect;
			bool intersectedFace = false;
			bool intersectedEdge = false;

			math::Point p_origin = p.distance(p_int) < p.distance(p_axe) ? p_int : p_axe;

			double dist_p = p.distance(p_axe);
			double dist_pm = p_m.distance(p_axe);

			double dist1;
			double dist2;
			double omega;


			if(dist_p < dist_pm){
				std::vector<Face> nearest_locked_faces;
				for(const auto& f : nearest_morphed.get<Face>()){
					if (f.getIDs<Region>().size() == 2) {
						if ((m_mesh->isMarked<Region>(f.getIDs<Region>()[0],m_morphRegions) && !m_mesh->isMarked<Region>(f.getIDs<Region>()[1],m_morphRegions))
						    || (!m_mesh->isMarked<Region>(f.getIDs<Region>()[0],m_morphRegions) && m_mesh->isMarked<Region>(f.getIDs<Region>()[1],m_morphRegions))){
							nearest_locked_faces.push_back(f);
						}
					}
				}
				for(const auto& f : nearest_locked_faces) {
					math::Triangle t1(f.get<Node>()[0].point(),f.get<Node>()[1].point(),f.get<Node>()[2].point());
					math::Triangle t2(f.get<Node>()[0].point(),f.get<Node>()[2].point(),f.get<Node>()[3].point());

					if(ray.intersect3D(t1, p_intersect)){
						nearest_locked_face = f;
						intersectedFace = true;
						break;
					}else if(ray.intersect3D(t2, p_intersect)) {
						nearest_locked_face = f;
						intersectedFace = true;
						break;
					}
				}
				if(!intersectedFace){
					for (const auto& e : nearest_morphed.get<Edge>()) {
						double paramSeg;
						double paramRay;
						if(ray.intersect3D(e.segment(),p_intersect, paramSeg,paramRay)){
							intersectedEdge = true;
							nearest_locked_edge = e;
							break;
						}
					}
				}
				if(!intersectedFace && !intersectedEdge){
					std::cout<<"ERROR"<<std::endl;
				}

				p_m = p_intersect;

				dist1 = p.distance(p_origin);
				dist2 = p_m.distance(p_origin);
				omega = dist1 / dist2;
			}else if(dist_p > dist_pm){

				for(const auto& f : outerSkin) {
					math::Triangle t1(f.get<Node>()[0].point(),f.get<Node>()[1].point(),f.get<Node>()[2].point());
					math::Triangle t2(f.get<Node>()[0].point(),f.get<Node>()[2].point(),f.get<Node>()[3].point());

					if(ray.intersect3D(t1, p_intersect)){
						nearest_locked_face = f;
						intersectedFace = true;
						break;
					}else if(ray.intersect3D(t2, p_intersect)) {
						nearest_locked_face = f;
						intersectedFace = true;
						break;
					}
				}
				if(!intersectedFace){
					for (const auto& e : nearest_morphed.get<Edge>()) {
						double paramSeg;
						double paramRay;
						if(ray.intersect3D(e.segment(),p_intersect, paramSeg,paramRay)){
							intersectedEdge = true;
							nearest_locked_edge = e;
							break;
						}
					}
				}
				if(!intersectedFace && !intersectedEdge){
					std::cout<<"ERROR"<<std::endl;
				}

				p_m = p_intersect;

				dist1 = p.distance(p_ext);
				dist2 = p_m.distance(p_ext);
				omega = dist1 / dist2;
			}else{
				std::cout<<"error where is p ?"<<std::endl;
			}

			math::Vector3d vec;
			vec.setXYZ(0,0,0);
			if(intersectedFace){

				std::vector<math::Point> points012 = {nearest_locked_face.get<Node>()[0].point(),nearest_locked_face.get<Node>()[1].point(),
				                                      nearest_locked_face.get<Node>()[2].point()};
				std::vector<math::Point> points023 = {nearest_locked_face.get<Node>()[0].point(),nearest_locked_face.get<Node>()[2].point(),
				                                      nearest_locked_face.get<Node>()[3].point()};
				std::vector<math::Point> points013 = {nearest_locked_face.get<Node>()[0].point(),nearest_locked_face.get<Node>()[1].point(),
				                                      nearest_locked_face.get<Node>()[3].point()};
				std::vector<math::Point> points123 = {nearest_locked_face.get<Node>()[1].point(),nearest_locked_face.get<Node>()[2].point(),
				                                      nearest_locked_face.get<Node>()[3].point()};


				math::Triangle t0(points012[0],points012[1],points012[2]);
				math::Triangle t1(points023[0],points023[1],points023[2]);
				math::Triangle t2(points013[0],points013[1],points013[2]);
				math::Triangle t3(points123[0],points123[1],points123[2]);

				std::vector<double> coefs012 = {0,0,0};

				math::Vector3d vec0;
				math::Vector3d vec1;
				math::Vector3d vec2;
				math::Vector3d vec3;

				math::Point refined_intersect;
				int nb_intersection = 0;
				if(ray.intersect3D(t0,refined_intersect)){
					math::Point::computeBarycentric(points012,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*(coefs012[2]);

					vec += vec0+vec1+vec2;
					nb_intersection++;
				}
				if(ray.intersect3D(t1,refined_intersect)){
					math::Point::computeBarycentric(points023,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*(coefs012[2]);
					vec += vec0+vec1+vec2;
					nb_intersection++;
				}
				if(ray.intersect3D(t2,refined_intersect)){
					math::Point::computeBarycentric(points013,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*(coefs012[2]);
					vec += vec0+vec1+vec2;
					nb_intersection++;
				}
				if(ray.intersect3D(t3,refined_intersect)){
					math::Point::computeBarycentric(points123,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*(coefs012[2]);
					vec += vec0+vec1+vec2;
					nb_intersection++;
				}

				vec = vec/nb_intersection;

				/*math::Vector3d vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*((coefs012[0]+coefs023[0])/2);
				math::Vector3d vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*coefs012[1];
				math::Vector3d vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*((coefs012[2]+coefs023[2])/2);
				math::Vector3d vec3 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*coefs023[2];*/

				//vec = m_vecs[nearest_morphed.id()];
				//vec = vec0+vec1+vec2;

			}else if(intersectedEdge){

				math::Vector3d vec0 = m_vecs[nearest_locked_edge.getIDs<Node>()[0]];
				math::Vector3d vec1 = m_vecs[nearest_locked_edge.getIDs<Node>()[1]];

				//vec = m_vecs[nearest_morphed.id()];
				vec = (vec0+vec1)/2;

			}
			vec *= omega;

			m_vecs[n.id()] = vec;
		}
	}

	std::vector<TCellID> node_list(internal_modified_nodes.begin(),internal_modified_nodes.end());

	return node_list;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> EllipticMorph::noExteriorLock(){

	FastLocalize fl_int(m_lockedIn);
	FastLocalize fl_morphed(m_morphed);

	std::set<TCellID> internal_modified_nodes;

	for(int i = 0; i < m_ellipses.size()-1; i++) {

		double Xmin = m_ellipses[i][0];
		double Xmax = m_ellipses[i + 1][0];

		std::vector<Node> interval_morphed;

		for (auto n : m_morphed) {
			if (n.X() >= Xmin && n.X() < Xmax) {
				interval_morphed.push_back(n);
			}
		}

		double coefa_min = m_ellipses[i][1];
		double coefb_min = m_ellipses[i][2];
		double coefc_min = m_ellipses[i][3];

		double coefa_max = m_ellipses[i + 1][1];
		double coefb_max = m_ellipses[i + 1][2];
		double coefc_max = m_ellipses[i + 1][3];

		for (auto n : interval_morphed) {
			math::Point p = n.point();

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();

			double distXmin = p.X() - Xmin;
			double distXmax = Xmax - Xmin;
			double w = distXmin / distXmax;

			double coefy_min = p.Y() >= 0 ? coefa_min : coefc_min;
			double coefy_max = p.Y() >= 0 ? coefa_max : coefc_max;

			double coefy_test = coefy_min * (3 * pow(1 - w, 2) - 2 * pow(1 - w, 3)) + coefy_max * (3 * pow(w, 2) - 2 * pow(w, 3));
			double coefz_test = coefb_min * (3 * pow(1 - w, 2) - 2 * pow(1 - w, 3)) + coefb_max * (3 * pow(w, 2) - 2 * pow(w, 3));

			double theta = atan(n.Y() / n.Z());
			if (theta < 0) {
				theta *= -1;
			}

			double sinus = sin(theta);
			double cosinus = cos(theta);

			double newY = (p.Y() * ((coefy_test * sinus) + coefz_test * (1 - sinus)));
			double newZ = (p.Z() * ((coefz_test * cosinus) + coefy_test * (1 - cosinus)));

			math::Vector3d vec;
			vec.setXYZ(0, newY - p.Y(), newZ - p.Z());
			m_vecs[n.id()] = vec;
			internal_modified_nodes.insert(n.id());
		}
	}

	for(int i = 0; i < m_ellipses.size()-1; i++) {

		double Xmin = m_ellipses[i][0];
		double Xmax = m_ellipses[i+1][0];

		std::vector<Node> execution;

		for(auto n : m_mesh->nodes()){

			Node node = m_mesh->get<Node>(n);

			if (node.X() >= Xmin && node.X() < Xmax && !m_mesh->isMarked<Node>(n, m_locked)) {
				execution.push_back(node);
			}
		}

		for (const auto& n : execution) {
			internal_modified_nodes.insert(n.id());

			math::Point p = n.point();

			Node nearest_morphed = m_mesh->get<Node>(fl_morphed.find(p));

			std::vector<Face> nearest_locked_faces;
			for(const auto& f : nearest_morphed.get<Face>()){
				if (f.getIDs<Region>().size() == 2) {
					if ((m_mesh->isMarked<Region>(f.getIDs<Region>()[0],m_morphRegions) && !m_mesh->isMarked<Region>(f.getIDs<Region>()[1],m_morphRegions))
					    || (!m_mesh->isMarked<Region>(f.getIDs<Region>()[0],m_morphRegions) && m_mesh->isMarked<Region>(f.getIDs<Region>()[1],m_morphRegions))){
						nearest_locked_faces.push_back(f);
					}
				}
			}

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			math::Point p_m = nearest_morphed.point();
			math::Point p_axe(p.X(), 0, 0);

			math::Ray ray(p_axe, p);

			Face nearest_locked_face;
			Edge nearest_locked_edge;
			math::Point p_intersect;
			bool intersectedFace = false;
			bool intersectedEdge = false;
			for(const auto& f : nearest_locked_faces) {
				math::Triangle t1(f.get<Node>()[0].point(),f.get<Node>()[1].point(),f.get<Node>()[2].point());
				math::Triangle t2(f.get<Node>()[0].point(),f.get<Node>()[2].point(),f.get<Node>()[3].point());

				if(ray.intersect3D(t1, p_intersect)){
					nearest_locked_face = f;
					intersectedFace = true;
					break;
				}else if(ray.intersect3D(t2, p_intersect)) {
					nearest_locked_face = f;
					intersectedFace = true;
					break;
				}
			}
			if(!intersectedFace){
				for (const auto& e : nearest_morphed.get<Edge>()) {
					double paramSeg;
					double paramRay;
					if(ray.intersect3D(e.segment(),p_intersect, paramSeg,paramRay)){
						intersectedEdge = true;
						nearest_locked_edge = e;
						break;
					}
				}
			}
			if(!intersectedFace && !intersectedEdge){
				std::cout<<"ERROR"<<std::endl;
			}

			p_m = p_intersect;

			math::Point p_origin = p.distance(p_int) < p.distance(p_axe) ? p_int : p_axe;


			double dist_p = p.distance(p_axe);
			double dist_pm = p_m.distance(p_axe);

			double dist1;
			double dist2;
			double omega;


			if(dist_p < dist_pm){
				dist1 = p.distance(p_origin);
				dist2 = p_m.distance(p_origin);
				omega = dist1 / dist2;
			}else{
				std::cout<<"error where is p ?"<<std::endl;
			}

			//std::vector<TCellID> test = nearest_locked_face.getIDs<Node>();

			math::Vector3d vec;
			vec.setXYZ(0,0,0);
			if(intersectedFace){

				std::vector<math::Point> points012 = {nearest_locked_face.get<Node>()[0].point(),nearest_locked_face.get<Node>()[1].point(),
				   											nearest_locked_face.get<Node>()[2].point()};
				std::vector<math::Point> points023 = {nearest_locked_face.get<Node>()[0].point(),nearest_locked_face.get<Node>()[2].point(),
				                                      nearest_locked_face.get<Node>()[3].point()};
				std::vector<math::Point> points013 = {nearest_locked_face.get<Node>()[0].point(),nearest_locked_face.get<Node>()[1].point(),
				                                      nearest_locked_face.get<Node>()[3].point()};
				std::vector<math::Point> points123 = {nearest_locked_face.get<Node>()[1].point(),nearest_locked_face.get<Node>()[2].point(),
				                                      nearest_locked_face.get<Node>()[3].point()};


				math::Triangle t0(points012[0],points012[1],points012[2]);
				math::Triangle t1(points023[0],points023[1],points023[2]);
				math::Triangle t2(points013[0],points013[1],points013[2]);
				math::Triangle t3(points123[0],points123[1],points123[2]);

				std::vector<double> coefs012 = {0,0,0};

				math::Vector3d vec0;
				math::Vector3d vec1;
				math::Vector3d vec2;
				math::Vector3d vec3;

				math::Point refined_intersect;
				int nb_intersection = 0;
				if(ray.intersect3D(t0,refined_intersect)){
					math::Point::computeBarycentric(points012,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*(coefs012[2]);

					vec += vec0+vec1+vec2;
					nb_intersection++;
				}
				if(ray.intersect3D(t1,refined_intersect)){
					math::Point::computeBarycentric(points023,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*(coefs012[2]);
					vec += vec0+vec1+vec2;
					nb_intersection++;
				}
				if(ray.intersect3D(t2,refined_intersect)){
					math::Point::computeBarycentric(points013,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*(coefs012[2]);
					vec += vec0+vec1+vec2;
					nb_intersection++;
				}
				if(ray.intersect3D(t3,refined_intersect)){
					math::Point::computeBarycentric(points123,p_intersect,coefs012);

					vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*(coefs012[0]);
					vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*coefs012[1];
					vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*(coefs012[2]);
					vec += vec0+vec1+vec2;
					nb_intersection++;
				}

				vec = vec/nb_intersection;

				/*math::Vector3d vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*((coefs012[0]+coefs023[0])/2);
				math::Vector3d vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*coefs012[1];
				math::Vector3d vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*((coefs012[2]+coefs023[2])/2);
				math::Vector3d vec3 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*coefs023[2];*/

				//vec = m_vecs[nearest_morphed.id()];
				//vec = vec0+vec1+vec2;

			}else if(intersectedEdge){

				math::Vector3d vec0 = m_vecs[nearest_locked_edge.getIDs<Node>()[0]];
				math::Vector3d vec1 = m_vecs[nearest_locked_edge.getIDs<Node>()[1]];

				//vec = m_vecs[nearest_morphed.id()];
				vec = (vec0+vec1)/2;

			}
			vec *= omega;


			/*if(n.id() == 22045 || n.id() == 22054 || n.id() == 22063 || n.id() == 74843 || n.id() == 74851 || n.id() == 74859 || n.id() == 74844
			    || n.id() == 74852 || n.id() == 74860){
				std::cout<<"node "<<n<<" -> "<<n.point()+vec<<", (/\\) = "<<omega<<std::endl;
				std::cout<<"locked "<<nearest_morphed<<std::endl;
				if(n.id() == 74851 || n.id() == 74852){
					std::cout<<"origin vec "<<n.point()<<" -> "<<n.point() + m_vecs[nearest_morphed.id()]<<std::endl;
				}
			}*/

			m_vecs[n.id()] = vec;
		}
	}

	std::vector<TCellID> node_list(internal_modified_nodes.begin(),internal_modified_nodes.end());

	return node_list;
}
/*----------------------------------------------------------------------------*/
void EllipticMorph::markLockedCells()
{
	std::set<TCellID> set_li;
	std::set<TCellID> set_lo;
	std::set<TCellID> set_mo;

	var_morhR = m_mesh->newVariable<int,GMDS_REGION>("in_morph");

	if (m_inputName.find(".vtk") != -1) {
		//std::vector<Variable<int>* > morphed_var_list;
		for(auto const& name : m_to_morph){
			Variable<int>* current_morph_var = m_mesh->getVariable<int, GMDS_REGION>(name);
			for(auto r : m_mesh->regions()){
				Region region = m_mesh->get<Region>(r);
				if(current_morph_var->value(r) == 1)
					for(auto n : region.getIDs<Node>()) {
						m_mesh->mark<Region>(r,m_morphRegions);
						var_morhR->set(r, 1);
						set_mo.insert(n);
					}
			}
		}
		for(auto const& name : m_ext_lock){
			Variable<int>* current_morph_var = m_mesh->getVariable<int, GMDS_REGION>(name);
			for(auto r : m_mesh->regions()){
				Region region = m_mesh->get<Region>(r);
				if(current_morph_var->value(r) == 1)
					for(auto n : region.getIDs<Node>())
						set_lo.insert(n);
			}
		}
		for(auto const& name : m_int_lock){
			Variable<int>* current_morph_var = m_mesh->getVariable<int, GMDS_REGION>(name);
			for(auto r : m_mesh->regions()){
				Region region = m_mesh->get<Region>(r);
				if(current_morph_var->value(r) == 1)
					for(auto n : region.getIDs<Node>())
						set_li.insert(n);
			}
		}
	}
	else if(m_inputName.find(".mli") != -1){
		for(auto const& name : m_to_morph){
			CellGroup<Region>* group = m_mesh->getGroup<Region>(name);
			for(auto r : group->cells()){
				Region region = m_mesh->get<Region>(r);
				for(auto n : region.getIDs<Node>()) {
					m_mesh->mark<Region>(r, m_morphRegions);
					set_mo.insert(n);
				}
			}
		}
		for(auto const& name : m_ext_lock){
			CellGroup<Region>* group = m_mesh->getGroup<Region>(name);
			for(auto r : group->cells()){
				Region region = m_mesh->get<Region>(r);
				for(auto n : region.getIDs<Node>())
					set_lo.insert(n);
			}
		}
		for(auto const& name : m_int_lock){
			CellGroup<Region>* group = m_mesh->getGroup<Region>(name);
			for(auto r : group->cells()){
				Region region = m_mesh->get<Region>(r);
				for(auto n : region.getIDs<Node>())
					set_li.insert(n);
			}
		}
	}
	else{
		throw (GMDSException("Error file format not supported"));
	}

	for (auto n : set_li) {
		Node node = m_mesh->get<Node>(n);
		m_lockedIn.push_back(node);
		m_mesh->mark<Node>(n, m_locked);
	}
	for (auto n : set_lo) {
		Node node = m_mesh->get<Node>(n);
		m_lockedOut.push_back(node);
		m_mesh->mark<Node>(n, m_locked);
	}
	for (auto n : set_mo) {
		Node node = m_mesh->get<Node>(n);
		m_morphed.push_back(node);
		m_mesh->mark<Node>(n, m_locked);
	}
	for (auto n : m_mesh->nodes()) {
		Node node = m_mesh->get<Node>(n);
		if (node.Y() < 0.1 && node.Y() > -0.1 && node.Z() < 0.1 && node.Z() > -0.1) {
			m_mesh->mark(node, m_locked);
		}
	}

	Variable<int>* var_int = m_mesh->getVariable<int, GMDS_REGION>("int");
	Variable<int>* var_couche = m_mesh->getVariable<int, GMDS_REGION>("couche");

	/*for(auto f : m_mesh->faces()) {
		Face face = m_mesh->get<Face>(f);
		if (face.getIDs<Region>().size() == 2) {
			if (var_int->value(face.getIDs<Region>()[0]) == 1 && var_couche->value(face.getIDs<Region>()[1]) == 1
			    || var_int->value(face.getIDs<Region>()[1]) == 1 && var_couche->value(face.getIDs<Region>()[0]) == 1){
				m_mesh->mark<Face>(f, m_innerSkin);
			}
		}
	}*/


}
/*----------------------------------------------------------------------------*/
void EllipticMorph::finalize(){

	std::cout<<"------------- Writing result -------------"<<std::endl;

	if(m_outputName.find(".vtk") != -1) {
		std::cout<<"\t -> Write VTK format"<<std::endl;
		gmds::IGMeshIOService ioService(m_mesh);
		gmds::VTKWriter w(&ioService);
		w.setCellOptions(gmds::N | gmds::R);
		w.setDataOptions(gmds::N | gmds::R);
		w.write(m_outputName);
	}else if(m_inputName.find(".mli") != -1) {
#ifdef WITH_LIMA
		std::cout<<"\t -> Write Lima format"<<std::endl;
		gmds::LimaWriter limaWriter(*m_mesh);
		limaWriter.write(m_outputName, R|N);
#endif     // WITH_LIMA
	}else{
		throw GMDSException("File format not supported, only .vtk or .mli2 supported");
	}
}
