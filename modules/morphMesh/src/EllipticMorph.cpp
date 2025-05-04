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
#include <limits>

#ifdef WITH_LIMA
	#include "gmds/io/LimaReader.h"
	#include "gmds/io/LimaWriter.h"
#endif //WITH_LIMA
#	include "gmds/math/Ray.h"
#	include "gmds/math/Triangle.h"

/*----------------------------------------------------------------------------*/

using namespace gmds;
using namespace morphmesh;

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

		   double X, coefa, coefb, coefc, coefd, coefe;

		   X = std::stod(current_word);
		   *m_stream >> current_word;
		   coefa = std::stod(current_word);
		   *m_stream >> current_word;
		   coefb = std::stod(current_word);
		   *m_stream >> current_word;
		   coefc = std::stod(current_word);
		   *m_stream >> current_word;
		   coefd = std::stod(current_word);
		   *m_stream >> current_word;
		   coefe = std::stod(current_word);

		   m_ellipses.push_back({X,coefa,coefb,coefc,coefd,coefe});

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

	   m_minEdgeLength = std::numeric_limits<double>::max();
	   for(const auto& e : m_mesh->edges()){
		   Edge edge = m_mesh->get<Edge>(e);
		   if(edge.length() < m_minEdgeLength){
			   m_minEdgeLength = edge.length();
		   }
	   }

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

	std::cout<<" ------------------ Computing translations ended ------------------ "<<std::endl;

	for(auto n : modified_nodes){

		math::Vector3d vec(m_vecs[n]);
		Node node = m_mesh->get<Node>(n);
		double y = node.Y();
		double z = node.Z();

		node.setY(y+vec.Y());
		node.setZ(z+vec.Z());
	}

   std::cout<<" ------------------ Correct nodes translated ------------------ "<<std::endl;

	while(!m_not_treated_nodes.empty()){
		std::vector<TCellID> next_iteration;
      for(auto n_id : m_not_treated_nodes) {
			Node n = m_mesh->get<Node>(n_id);
			std::set<TCellID> neighbor_nodes;
			for (const auto &r : n.get<Region>()) {
				for (const auto &n_r : r.getIDs<Node>()) {
					if (n_r != n_id) neighbor_nodes.insert(n_r);
				}
			}
			math::Vector3d vec;
			vec.setXYZ(0, 0, 0);
			int cpt = 0;
			for (auto node : neighbor_nodes) {

				if (m_vecs.find(node) != m_vecs.end()) {
					vec = vec + m_vecs[node];
					cpt++;
				}
			}
			if(cpt > 0 ){
            vec = vec / cpt;

            m_vecs[n.id()] = vec;

            double y = n.Y();
            double z = n.Z();


            n.setY(y + vec.Y());
            n.setZ(z + vec.Z());
			}
			else{
				next_iteration.push_back(n_id);
			}
		}

		m_not_treated_nodes = next_iteration;
	}

   std::cout<<" ------------------ Incorrect nodes translated ------------------ "<<std::endl;


   finalize();

	std::cout<<"============================== Elliptic Morphing Finished =============================="<<std::endl;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> EllipticMorph::withExteriorLock(){

	Variable<int>* treated = m_mesh->getVariable<int, GMDS_NODE>("not_treated");

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
			if (n.X() >= Xmin && n.X() < Xmax && n.Z() >=0) {
				interval_morphed.push_back(n);
			}
		}

		double coefa_min = m_ellipses[i][1];
		double coefb_min = m_ellipses[i][2];
		double coefc_min = m_ellipses[i][3];
		double coefd_min = m_ellipses[i][4];
		double coefe_min = m_ellipses[i][5];

		double coefa_max = m_ellipses[i + 1][1];
		double coefb_max = m_ellipses[i + 1][2];
		double coefc_max = m_ellipses[i + 1][3];
		double coefd_max = m_ellipses[i + 1][4];
		double coefe_max = m_ellipses[i + 1][5];

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
			double coef45_test = coefd_min * (3 * pow(1 - w, 2) - 2 * pow(1 - w, 3)) + coefd_max * (3 * pow(w, 2) - 2 * pow(w, 3));
			double coef135_test = coefe_min * (3 * pow(1 - w, 2) - 2 * pow(1 - w, 3)) + coefe_max * (3 * pow(w, 2) - 2 * pow(w, 3));

			math::Vector3d vY;
			vY.setXYZ(0,1,0);

			math::Vector3d coords;
			coords.setXYZ(0,p.Y(),p.Z());

			double dot = vY.dot(coords);

			double theta = atan(n.Y()/n.Z());
			if (theta < 0) {
				//theta *= -1;
			}
			//std::cout<<theta<<std::endl;

			double sinus = sin(theta);
			double cosinus = cos(theta);

			if (theta < 0) {
				sinus *= -1;
				//cosinus *= -1;
			}

			double coef1, coef2;
			if(theta <= math::Constants::PIDIV2 && theta >= math::Constants::PIDIV4){
				coef1 = coefy_test;
				coef2 = coef45_test;

				sinus = (sinus - (sqrt(2)/2))/ (1-(sqrt(2)/2));
				cosinus = cosinus / (sqrt(2)/2);
			}else if(theta < math::Constants::PIDIV4 && theta >= 0){
				coef1 = coef45_test;
				coef2 = coefz_test;

				sinus = sinus / (sqrt(2)/2);
				cosinus = (cosinus - (sqrt(2)/2))/ (1-(sqrt(2)/2));
			}else if(theta > -math::Constants::PIDIV4 && theta <0){
				coef1 = coef135_test;
				coef2 = coefz_test;

				sinus = sinus / (sqrt(2)/2);
				cosinus = (cosinus - (sqrt(2)/2))/ (1-(sqrt(2)/2));
			}else if(theta <= -math::Constants::PIDIV4 && theta >= -math::Constants::PIDIV2){
				coef1 = coefy_test;
				coef2 = coef135_test;

				sinus = (sinus - (sqrt(2)/2))/ (1-(sqrt(2)/2));
				cosinus = cosinus / (sqrt(2)/2);
			}else{
				std::cout<<"ERROR"<<std::endl;
			}


			double newY = (p.Y() * ((coef1 * sinus) + coef2 * (1 - sinus)));
			double newZ = (p.Z() * ((coef2 * cosinus) + coef1 * (1 - cosinus)));

			math::Vector3d vec;
			vec.setXYZ(0, newY - p.Y(), newZ - p.Z());
			m_vecs[n.id()] = vec;
			internal_modified_nodes.insert(n.id());
		}
	}

	for(int i = 0; i < m_ellipses.size()-1; i++) {

		double Xmin = m_ellipses[i][0];
		double Xmax = m_ellipses[i+1][0];

		std::set<Face> faces_morphed_in_ellipse_zone;
		for(auto r : m_mesh->regions()) {
			if (m_mesh->isMarked<Region>(r, m_morphRegions)) {
				std::vector<Face> f_r = m_mesh->get<Region>(r).get<Face>();
				for (const auto &f : f_r) {

					int cptIsIn = 4;
					for (const auto &n_f : f.get<Node>()) {
						if (n_f.X() < Xmin - m_minEdgeLength || n_f.X() > Xmax + m_minEdgeLength) {
							cptIsIn--;
						}
					}

					if (cptIsIn > 0) {
						if (f.getIDs<Region>().size() == 2) {
							if ((m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions) && !m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))
							    || (!m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions)
							        && m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))) {

								faces_morphed_in_ellipse_zone.insert(f);
							}
						}
					}
				}
			}
		}

		std::set<Face> faces_lock_in_ellipse_zone;
		for(auto f : outerSkin) {

			int cptIsIn = 4;
			for (const auto &n_f : f.get<Node>()) {
				if (n_f.X() < Xmin - m_minEdgeLength || n_f.X() > Xmax + m_minEdgeLength) {
					cptIsIn--;
				}
			}

			if (cptIsIn > 0) {
				faces_lock_in_ellipse_zone.insert(f);
			}
		}

		std::vector<Node> execution;

		for(auto n : m_mesh->nodes()){

			Node node = m_mesh->get<Node>(n);

			if (node.X() >= Xmin && node.X() < Xmax && !m_mesh->isMarked<Node>(n, m_locked)  && node.Z() >=0) {
				execution.push_back(node);
			}
		}

		for (const auto& n : execution) {
			math::Point p = n.point();

			Node nearest_morphed = m_mesh->get<Node>(fl_morphed.find(p));

			// les faces adjacentes au morphed point le plus proche pour réduire les faces à intersecter
			std::vector<Face> nearest_morphed_faces;
			for (const auto &f : nearest_morphed.get<Face>()) {
				if (f.getIDs<Region>().size() == 2) {
					if ((m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions) && !m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))
					    || (!m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions) && m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))) {
						nearest_morphed_faces.push_back(f);
					}
				}
			}
			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			math::Point p_m = nearest_morphed.point();
			math::Point p_axe(p.X(), 0, 0);
			math::Point p_ext;

			math::Ray ray(p_axe, p);


			Face nearest_morphed_face;
			Edge nearest_morphed_edge;
			math::Point p_intersect;
			bool intersectedFace = false;
			bool intersectedEdge = false;

			//-------------------------------------------------------------------------
			// Step to identified intersection point on the morphed surface
			for (const auto &f : nearest_morphed_faces) {
				math::Triangle t1(f.get<Node>()[0].point(), f.get<Node>()[1].point(), f.get<Node>()[2].point());
				math::Triangle t2(f.get<Node>()[0].point(), f.get<Node>()[2].point(), f.get<Node>()[3].point());

				if (ray.intersect3D(t1, p_intersect)) {
					nearest_morphed_face = f;
					intersectedFace = true;
					break;
				}
				else if (ray.intersect3D(t2, p_intersect)) {
					nearest_morphed_face = f;
					intersectedFace = true;
					break;
				}
			}
			if (!intersectedFace) {

				std::set<TCellID> neighbors;
				for (const auto &f : nearest_morphed_faces) {
					for (auto n_f : f.getIDs<Node>()) {
						neighbors.insert(n_f);
					}
				}
				for (const auto &e : nearest_morphed.get<Edge>()) {
					bool correct_edge = false;
					if (e.getIDs<Node>()[0] == nearest_morphed.id()) {
						if (neighbors.find(e.getIDs<Node>()[1]) != neighbors.end()) {
							correct_edge = true;
						}
					}
					else if (e.getIDs<Node>()[1] == nearest_morphed.id()) {
						if (neighbors.find(e.getIDs<Node>()[0]) != neighbors.end()) {
							correct_edge = true;
						}
					}
					if (correct_edge) {
						double paramSeg;
						double paramRay;
						if (ray.intersect3D(e.segment(), p_intersect, paramSeg, paramRay)) {
							intersectedEdge = true;
							nearest_morphed_edge = e;
							break;
						}
					}
				}
			}
			bool noIntersection = false;
			if (!intersectedFace && !intersectedEdge) {
				noIntersection = true;
				// std::cout<<"ERROR"<<std::endl;

				std::vector<Face> intersectedFaces;
				std::vector<math::Point> points;
				for (const auto &f : faces_morphed_in_ellipse_zone) {
					int cptIsIn = 4; //On prend que les faces qui sont dans la meme "tranche" X que le noeud n pour pas tout tester
					for (const auto &n_f : f.get<Node>()) {
						if (n_f.X() < n.X() - m_minEdgeLength || n_f.X() > n.X() + m_minEdgeLength) {
							cptIsIn--;
						}
					}

					if (cptIsIn > 0) {
						math::Triangle t1(f.get<Node>()[0].point(), f.get<Node>()[1].point(), f.get<Node>()[2].point());
						math::Triangle t2(f.get<Node>()[0].point(), f.get<Node>()[2].point(), f.get<Node>()[3].point());

						if (ray.intersect3D(t1, p_intersect)) {
							intersectedFaces.push_back(f);
							points.push_back(p_intersect);
							intersectedFace = true;
							break;
						}
						else if (ray.intersect3D(t2, p_intersect)) {
							intersectedFaces.push_back(f);
							points.push_back(p_intersect);
							intersectedFace = true;
							break;
						}
					}

					if (!intersectedFaces.empty()) {
						double distance = std::numeric_limits<double>::max();
						for (int i_f = 0; i_f < intersectedFaces.size(); i_f++) {
							if (points[i_f].distance(p) < distance) {
								distance = points[i_f].distance(p);
								nearest_morphed_face = intersectedFaces[i_f];
								p_intersect = points[i_f];
							}
						}
						noIntersection = false;
					}
					// We may need to add the case where the ray intersect an Edge
				}
				if(!intersectedFace) {
					for (const auto &f : faces_lock_in_ellipse_zone) {
						int cptIsIn = 4;     // On prend que les faces qui sont dans la meme "tranche" X que le noeud n pour pas tout tester
						for (const auto &n_f : f.get<Node>()) {
							if (n_f.X() < n.X() - m_minEdgeLength || n_f.X() > n.X() + m_minEdgeLength) {
								cptIsIn--;
							}
						}

						if (cptIsIn > 0) {

							std::vector<Edge> edges = f.get<Edge>();
							for(auto e : edges) {
								double paramSeg;
								double paramRay;
								if (ray.intersect3D(e.segment(), p_intersect, paramSeg, paramRay)) {
									intersectedEdge = true;
									nearest_morphed_edge = e;
									noIntersection = false;
									break;
								}
							}
						}
					}
				}
			}
			//-------------------------------------------------------------------------

			//-------------------------------------------------------------------------
			// Step to identified intersection point on the ext locked surface
			Face nearest_locked_face;
			Edge nearest_locked_edge;
			bool intersectedExtFace = false;
			bool intersectedExtEdge = false;
			bool noIntersectionExt = false;

			std::vector<Face> intersectedFaces;
			std::vector<math::Point> points;
			for (const auto &f : faces_lock_in_ellipse_zone) {
				int cptIsIn = 4; //On prend que les faces qui sont dans la meme "tranche" X que le noeud n pour pas tout tester
				for (const auto &n_f : f.get<Node>()) {
					if (n_f.X() < n.X() - m_minEdgeLength || n_f.X() > n.X() + m_minEdgeLength) {
						cptIsIn--;
					}
				}

				if (cptIsIn > 0) {
					math::Triangle t1(f.get<Node>()[0].point(), f.get<Node>()[1].point(), f.get<Node>()[2].point());
					math::Triangle t2(f.get<Node>()[0].point(), f.get<Node>()[2].point(), f.get<Node>()[3].point());

					if (ray.intersect3D(t1, p_ext)) {
						intersectedFaces.push_back(f);
						points.push_back(p_ext);
						intersectedExtFace = true;
						break;
					}
					else if (ray.intersect3D(t2, p_ext)) {
						intersectedFaces.push_back(f);
						points.push_back(p_ext);
						intersectedExtFace = true;
						break;
					}
				}

				if (!intersectedFaces.empty()) {
					double distance = std::numeric_limits<double>::max();
					for (int i_f = 0; i_f < intersectedFaces.size(); i_f++) {
						if (points[i_f].distance(p) < distance) {
							distance = points[i_f].distance(p);
							nearest_locked_face = intersectedFaces[i_f];
							p_ext = points[i_f];
						}
					}
					noIntersectionExt = false;
				}




				// We may need to add the case where the ray intersect an Edge
			}
		   if(!intersectedFace) {
				for (const auto &f : faces_lock_in_ellipse_zone) {
					int cptIsIn = 4;     // On prend que les faces qui sont dans la meme "tranche" X que le noeud n pour pas tout tester
					for (const auto &n_f : f.get<Node>()) {
						if (n_f.X() < n.X() - m_minEdgeLength || n_f.X() > n.X() + m_minEdgeLength) {
							cptIsIn--;
						}
					}

					if (cptIsIn > 0) {

						std::vector<Edge> edges = f.get<Edge>();
						for(auto e : edges) {
							double paramSeg;
							double paramRay;
							if (ray.intersect3D(e.segment(), p_intersect, paramSeg, paramRay)) {
								intersectedEdge = true;
								nearest_locked_edge = e;
								noIntersectionExt = false;
								break;
							}
						}
					}
				}
			}
			//-------------------------------------------------------------------------



			if (noIntersection || noIntersectionExt) {
				m_not_treated_nodes.push_back(n.id());
				treated->set(n.id(),1);

			}else {
				p_m = p_intersect;

				math::Point p_origin = p.distance(p_int) < p.distance(p_axe) ? p_int : p_axe;

				double dist_p = p.distance(p_axe);
				double dist_pm = p_m.distance(p_axe);

				double dist1;
				double dist2;
				double omega;

				if (dist_p < dist_pm) {
					dist1 = p.distance(p_origin);
					dist2 = p_m.distance(p_origin);
					omega = dist1 / dist2;
				}
				else if (dist_p > dist_pm) {
					dist1 = p.distance(p_ext);
					dist2 = p_m.distance(p_ext);
					omega = dist1 / dist2;
				}
				else {
					std::cout << "error where is p ?" << std::endl;
					std::cout<<"dist p = "<<dist_p<<" et dist pm ="<<dist_pm<<std::endl;
					std::cout<<n<<std::endl;
					std::cout<<nearest_morphed<<std::endl;
				}

				math::Vector3d vec;
				vec.setXYZ(0, 0, 0);
				if (intersectedFace) {

					std::vector<math::Point> points012 = {nearest_morphed_face.get<Node>()[0].point(), nearest_morphed_face.get<Node>()[1].point(),
					                                      nearest_morphed_face.get<Node>()[2].point()};
					std::vector<math::Point> points023 = {nearest_morphed_face.get<Node>()[0].point(), nearest_morphed_face.get<Node>()[2].point(),
					                                      nearest_morphed_face.get<Node>()[3].point()};
					std::vector<math::Point> points013 = {nearest_morphed_face.get<Node>()[0].point(), nearest_morphed_face.get<Node>()[1].point(),
					                                      nearest_morphed_face.get<Node>()[3].point()};
					std::vector<math::Point> points123 = {nearest_morphed_face.get<Node>()[1].point(), nearest_morphed_face.get<Node>()[2].point(),
					                                      nearest_morphed_face.get<Node>()[3].point()};

					math::Triangle t0(points012[0], points012[1], points012[2]);
					math::Triangle t1(points023[0], points023[1], points023[2]);
					math::Triangle t2(points013[0], points013[1], points013[2]);
					math::Triangle t3(points123[0], points123[1], points123[2]);

					std::vector<double> coefs012 = {0, 0, 0};

					math::Vector3d vec0;
					math::Vector3d vec1;
					math::Vector3d vec2;
					math::Vector3d vec3;

					math::Point refined_intersect;
					int nb_intersection = 0;
					if (ray.intersect3D(t0, refined_intersect)) {
						math::Point::computeBarycentric(points012, p_intersect, coefs012);

						vec0 = m_vecs[nearest_morphed_face.getIDs<Node>()[0]] * (coefs012[0]);
						vec1 = m_vecs[nearest_morphed_face.getIDs<Node>()[1]] * coefs012[1];
						vec2 = m_vecs[nearest_morphed_face.getIDs<Node>()[2]] * (coefs012[2]);

						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}
					if (ray.intersect3D(t1, refined_intersect)) {
						math::Point::computeBarycentric(points023, p_intersect, coefs012);

						vec0 = m_vecs[nearest_morphed_face.getIDs<Node>()[0]] * (coefs012[0]);
						vec1 = m_vecs[nearest_morphed_face.getIDs<Node>()[2]] * coefs012[1];
						vec2 = m_vecs[nearest_morphed_face.getIDs<Node>()[3]] * (coefs012[2]);
						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}
					if (ray.intersect3D(t2, refined_intersect)) {
						math::Point::computeBarycentric(points013, p_intersect, coefs012);

						vec0 = m_vecs[nearest_morphed_face.getIDs<Node>()[0]] * (coefs012[0]);
						vec1 = m_vecs[nearest_morphed_face.getIDs<Node>()[1]] * coefs012[1];
						vec2 = m_vecs[nearest_morphed_face.getIDs<Node>()[3]] * (coefs012[2]);
						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}
					if (ray.intersect3D(t3, refined_intersect)) {
						math::Point::computeBarycentric(points123, p_intersect, coefs012);

						vec0 = m_vecs[nearest_morphed_face.getIDs<Node>()[1]] * (coefs012[0]);
						vec1 = m_vecs[nearest_morphed_face.getIDs<Node>()[2]] * coefs012[1];
						vec2 = m_vecs[nearest_morphed_face.getIDs<Node>()[3]] * (coefs012[2]);
						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}

					vec = vec / nb_intersection;

					/*math::Vector3d vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]]*((coefs012[0]+coefs023[0])/2);
					math::Vector3d vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]]*coefs012[1];
					math::Vector3d vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[2]]*((coefs012[2]+coefs023[2])/2);
					math::Vector3d vec3 = m_vecs[nearest_locked_face.getIDs<Node>()[3]]*coefs023[2];*/

					// vec = m_vecs[nearest_morphed.id()];
					// vec = vec0+vec1+vec2;
				}
				else if (intersectedEdge) {

					math::Vector3d vec0 = m_vecs[nearest_morphed_edge.getIDs<Node>()[0]];
					math::Vector3d vec1 = m_vecs[nearest_morphed_edge.getIDs<Node>()[1]];

					// vec = m_vecs[nearest_morphed.id()];
					vec = (vec0 + vec1) / 2;
				}
				vec *= omega;
				m_vecs[n.id()] = vec;
				internal_modified_nodes.insert(n.id());
			}
		}
	}

	std::vector<TCellID> node_list(internal_modified_nodes.begin(),internal_modified_nodes.end());

	return node_list;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> EllipticMorph::noExteriorLock(){

	Variable<int>* treated = m_mesh->getVariable<int, GMDS_NODE>("not_treated");

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

      std::set<Face> faces_morphed_in_ellipse_zone;
      for(auto r : m_mesh->regions()) {
			if (m_mesh->isMarked<Region>(r, m_morphRegions)) {
				std::vector<Face> f_r = m_mesh->get<Region>(r).get<Face>();
				for (const auto &f : f_r) {

					int cptIsIn = 4;
					for (const auto &n_f : f.get<Node>()) {
						if (n_f.X() < Xmin - m_minEdgeLength || n_f.X() > Xmax + m_minEdgeLength) {
							cptIsIn--;
						}
					}

					if (cptIsIn > 0) {
						if (f.getIDs<Region>().size() == 2) {
							if ((m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions) && !m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))
							    || (!m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions)
							        && m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))) {

								faces_morphed_in_ellipse_zone.insert(f);
							}
						}
					}
				}
			}
		}

		std::vector<Node> execution;

		for(auto n : m_mesh->nodes()){

			Node node = m_mesh->get<Node>(n);

			if (node.X() >= Xmin && node.X() < Xmax && !m_mesh->isMarked<Node>(n, m_locked)) {
				execution.push_back(node);
			}
		}

		for (const auto& n : execution) {
			math::Point p = n.point();

			Node nearest_morphed = m_mesh->get<Node>(fl_morphed.find(p));

			std::vector<Face> nearest_morphed_faces;
			for (const auto &f : nearest_morphed.get<Face>()) {
				if (f.getIDs<Region>().size() == 2) {
					if ((m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions) && !m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))
					    || (!m_mesh->isMarked<Region>(f.getIDs<Region>()[0], m_morphRegions) && m_mesh->isMarked<Region>(f.getIDs<Region>()[1], m_morphRegions))) {
						nearest_morphed_faces.push_back(f);
					}
				}
			}

			math::Point p_int;
			if (!m_lockedIn.empty()) {
				p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			}
			else {
				p_int.setXYZ(p.X(), 0, 0);
			}

			math::Point p_m = nearest_morphed.point();
			math::Point p_axe(p.X(), 0, 0);

			math::Ray ray(p_axe, p);

			Face nearest_locked_face;
			Edge nearest_locked_edge;
			math::Point p_intersect;
			bool intersectedFace = false;
			bool intersectedEdge = false;

			for (const auto &f : nearest_morphed_faces) {
				math::Triangle t1(f.get<Node>()[0].point(), f.get<Node>()[1].point(), f.get<Node>()[2].point());
				math::Triangle t2(f.get<Node>()[0].point(), f.get<Node>()[2].point(), f.get<Node>()[3].point());

				if (ray.intersect3D(t1, p_intersect)) {
					nearest_locked_face = f;
					intersectedFace = true;
					break;
				}
				else if (ray.intersect3D(t2, p_intersect)) {
					nearest_locked_face = f;
					intersectedFace = true;
					break;
				}
			}
			if (!intersectedFace) {

				std::set<TCellID> neighbors;
				for (const auto &f : nearest_morphed_faces) {
					for (auto n_f : f.getIDs<Node>()) {
						neighbors.insert(n_f);
					}
				}
				for (const auto &e : nearest_morphed.get<Edge>()) {
					bool correct_edge = false;
					if (e.getIDs<Node>()[0] == nearest_morphed.id()) {
						if (neighbors.find(e.getIDs<Node>()[1]) != neighbors.end()) {
							correct_edge = true;
						}
					}
					else if (e.getIDs<Node>()[1] == nearest_morphed.id()) {
						if (neighbors.find(e.getIDs<Node>()[0]) != neighbors.end()) {
							correct_edge = true;
						}
					}
					if (correct_edge) {
						double paramSeg;
						double paramRay;
						if (ray.intersect3D(e.segment(), p_intersect, paramSeg, paramRay)) {
							intersectedEdge = true;
							nearest_locked_edge = e;
							break;
						}
					}
				}
			}
			bool noIntersection = false;

			if (!intersectedFace && !intersectedEdge) {
				noIntersection = true;
				// std::cout<<"ERROR"<<std::endl;

				std::vector<Face> intersectedFaces;
				std::vector<math::Point> points;
				for (const auto &f : faces_morphed_in_ellipse_zone) {
					int cptIsIn = 4; //On prend que les faces qui sont dans la meme "tranche" X que le noeud n pour pas tout tester
					for (const auto &n_f : f.get<Node>()) {
						if (n_f.X() < n.X() - m_minEdgeLength || n_f.X() > n.X() + m_minEdgeLength) {
							cptIsIn--;
						}
					}

					if (cptIsIn > 0) {
						math::Triangle t1(f.get<Node>()[0].point(), f.get<Node>()[1].point(), f.get<Node>()[2].point());
						math::Triangle t2(f.get<Node>()[0].point(), f.get<Node>()[2].point(), f.get<Node>()[3].point());

						if (ray.intersect3D(t1, p_intersect)) {
							intersectedFaces.push_back(f);
							points.push_back(p_intersect);
							intersectedFace = true;
							break;
						}
						else if (ray.intersect3D(t2, p_intersect)) {
							intersectedFaces.push_back(f);
							points.push_back(p_intersect);
							intersectedFace = true;
							break;
						}
					}

					if (!intersectedFaces.empty()) {
						double distance = std::numeric_limits<double>::max();
						for (int i_f = 0; i_f < intersectedFaces.size(); i_f++) {
							if (points[i_f].distance(p) < distance) {
								distance = points[i_f].distance(p);
								nearest_locked_face = intersectedFaces[i_f];
								p_intersect = points[i_f];
							}
						}
						noIntersection = false;
					}
					// We may need to add the case where the ray intersect an Edge
				}
			}

			if (noIntersection) {
				m_not_treated_nodes.push_back(n.id());
				treated->set(n.id(),1);
			}
			else {   //Ici on a trouvé une intersection on va faire le calcul du déplacement
				p_m = p_intersect;

				math::Point p_origin = p.distance(p_int) < p.distance(p_axe) ? p_int : p_axe;

				double dist_p = p.distance(p_axe);
				double dist_pm = p_m.distance(p_axe);

				double dist1;
				double dist2;
				double omega;

				if (dist_p < dist_pm) {
					dist1 = p.distance(p_origin);
					dist2 = p_m.distance(p_origin);
					omega = dist1 / dist2;
				}
				else {
					std::cout << "error where is p ?" << std::endl;
					std::cout<<"dist p = "<<dist_p<<" et dist pm ="<<dist_pm<<std::endl;
					std::cout<<n<<std::endl;
					std::cout<<nearest_morphed<<std::endl;
				}

				// std::vector<TCellID> test = nearest_locked_face.getIDs<Node>();

				math::Vector3d vec;
				vec.setXYZ(0, 0, 0);
				if (intersectedFace) {

					std::vector<math::Point> points012 = {nearest_locked_face.get<Node>()[0].point(), nearest_locked_face.get<Node>()[1].point(),
					                                      nearest_locked_face.get<Node>()[2].point()};
					std::vector<math::Point> points023 = {nearest_locked_face.get<Node>()[0].point(), nearest_locked_face.get<Node>()[2].point(),
					                                      nearest_locked_face.get<Node>()[3].point()};
					std::vector<math::Point> points013 = {nearest_locked_face.get<Node>()[0].point(), nearest_locked_face.get<Node>()[1].point(),
					                                      nearest_locked_face.get<Node>()[3].point()};
					std::vector<math::Point> points123 = {nearest_locked_face.get<Node>()[1].point(), nearest_locked_face.get<Node>()[2].point(),
					                                      nearest_locked_face.get<Node>()[3].point()};

					math::Triangle t0(points012[0], points012[1], points012[2]);
					math::Triangle t1(points023[0], points023[1], points023[2]);
					math::Triangle t2(points013[0], points013[1], points013[2]);
					math::Triangle t3(points123[0], points123[1], points123[2]);

					std::vector<double> coefs012 = {0, 0, 0};

					math::Vector3d vec0;
					math::Vector3d vec1;
					math::Vector3d vec2;
					math::Vector3d vec3;

					math::Point refined_intersect;
					int nb_intersection = 0;
					if (ray.intersect3D(t0, refined_intersect)) {
						math::Point::computeBarycentric(points012, p_intersect, coefs012);

						vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]] * (coefs012[0]);
						vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]] * coefs012[1];
						vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[2]] * (coefs012[2]);

						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}
					if (ray.intersect3D(t1, refined_intersect)) {
						math::Point::computeBarycentric(points023, p_intersect, coefs012);

						vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]] * (coefs012[0]);
						vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[2]] * coefs012[1];
						vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]] * (coefs012[2]);
						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}
					if (ray.intersect3D(t2, refined_intersect)) {
						math::Point::computeBarycentric(points013, p_intersect, coefs012);

						vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[0]] * (coefs012[0]);
						vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[1]] * coefs012[1];
						vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]] * (coefs012[2]);
						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}
					if (ray.intersect3D(t3, refined_intersect)) {
						math::Point::computeBarycentric(points123, p_intersect, coefs012);

						vec0 = m_vecs[nearest_locked_face.getIDs<Node>()[1]] * (coefs012[0]);
						vec1 = m_vecs[nearest_locked_face.getIDs<Node>()[2]] * coefs012[1];
						vec2 = m_vecs[nearest_locked_face.getIDs<Node>()[3]] * (coefs012[2]);
						vec += vec0 + vec1 + vec2;
						nb_intersection++;
					}

					vec = vec / nb_intersection;
				}
				else if (intersectedEdge) {
					math::Vector3d vec0 = m_vecs[nearest_locked_edge.getIDs<Node>()[0]];
					math::Vector3d vec1 = m_vecs[nearest_locked_edge.getIDs<Node>()[1]];

					vec = (vec0 + vec1) / 2;
				}
					vec *= omega;
					m_vecs[n.id()] = vec;
					internal_modified_nodes.insert(n.id());
			}
			//end iteration for node n
		}
		//end iteration for ellipse
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

	m_mesh->newVariable<int, GMDS_NODE>("not_treated");


	if (m_inputName.find(".vtk") != -1) {
		//std::vector<Variable<int>* > morphed_var_list;
		for(auto const& name : m_to_morph){
			Variable<int>* current_morph_var = m_mesh->getVariable<int, GMDS_REGION>(name);
			for(auto r : m_mesh->regions()){
				Region region = m_mesh->get<Region>(r);
				if(current_morph_var->value(r) == 1)
					for(auto n : region.getIDs<Node>()) {
						m_mesh->mark<Region>(r,m_morphRegions);
						set_mo.insert(n);
					}
			}
		}
		for(auto const& name : m_ext_lock){
			Variable<int>* current_morph_var = m_mesh->getVariable<int, GMDS_REGION>(name);
			for(auto r : m_mesh->regions()){
				Region region = m_mesh->get<Region>(r);
				if(current_morph_var->value(r) == 1){
					m_mesh->mark<Region>(r, m_lockRegions);
					for(auto n : region.getIDs<Node>())
						set_lo.insert(n);
				}
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
/*----------------------------------------------------------------------------*/
void EllipticMorph::writeMesh(std::string AFilename) const
{
	std::cout<<"------------- PRINT ------------- "<<AFilename<<std::endl;

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKWriter w(&ioService);
	w.setCellOptions(gmds::N | gmds::R);
	w.setDataOptions(gmds::N | gmds::R);
	w.write(AFilename);
}
/*----------------------------------------------------------------------------*/