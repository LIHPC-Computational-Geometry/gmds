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
	   doc.updateUpwardConnectivity();

	   m_locked = m_mesh->newMark<Node>();
}
/*----------------------------------------------------------------------------*/
EllipticMorph::~EllipticMorph() = default;
/*----------------------------------------------------------------------------*/
void EllipticMorph::execute()
{

	std::cout<<"============================== Starting Elliptic Morphing =============================="<<std::endl;

	markLockedCells();

	std::vector<Node> nodes_int;
	std::vector<Node> nodes_ext;

	FastLocalize fl_int(m_lockedIn);
	FastLocalize fl_ext(m_lockedOut);
	FastLocalize fl_morphed(m_morphed);

	std::map<TCellID, math::Vector3d> vecs;
	std::set<TCellID> modified_nodes;

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
			vecs[n.id()] = vec;
			modified_nodes.insert(n.id());
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
			modified_nodes.insert(n.id());

			math::Point p = n.point();

			Node nearest_morphed = m_mesh->get<Node>(fl_morphed.find(p));

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			math::Point p_ext = m_mesh->get<Node>(fl_ext.find(p)).point();
			math::Point p_m = nearest_morphed.point();
			math::Point p_axe(p.X(), 0, 0);

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
			}else if(dist_p > dist_pm){
				dist1 = p.distance(p_ext);
				dist2 = p_m.distance(p_ext);
				omega = dist1 / dist2;
			}else{
				std::cout<<"error where is p ?"<<std::endl;
			}

			math::Vector3d vec;
			vec = vecs[nearest_morphed.id()];
			vec *= omega;

			vecs[n.id()] = vec;
		}
	}

	for(auto n : modified_nodes){

		math::Vector3d vec(vecs[n]);
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
void EllipticMorph::markLockedCells()
{
	std::set<TCellID> set_li;
	std::set<TCellID> set_lo;
	std::set<TCellID> set_mo;

	if (m_inputName.find(".vtk") != -1) {
		//std::vector<Variable<int>* > morphed_var_list;
		for(auto const& name : m_to_morph){
			Variable<int>* current_morph_var = m_mesh->getVariable<int, GMDS_REGION>(name);
			for(auto r : m_mesh->regions()){
				Region region = m_mesh->get<Region>(r);
				if(current_morph_var->value(r) == 1)
					for(auto n : region.getIDs<Node>())
						set_mo.insert(n);
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
				for(auto n : region.getIDs<Node>())
					set_mo.insert(n);
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
