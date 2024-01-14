/*----------------------------------------------------------------------------*/
#include "gmds/morphMesh/MorphMesh.h"
#include "gmds/math/Plane.h"
#include "gmds/math/Triangle.h"
#include <iomanip>

/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace morphmesh {
/*----------------------------------------------------------------------------*/
MorphMesh::MorphMesh(Mesh *AMesh, const std::vector<math::Point> &APoints, double ARadius): m_mesh(AMesh)
{
	m_radius = ARadius;
	m_targets = APoints;
	m_locked = m_mesh->newMark<Node>();

	locked_faces = m_mesh->newVariable<int,GMDS_FACE>("locked_face");
}
/*----------------------------------------------------------------------------*/
MorphMesh::~MorphMesh() {}
/*----------------------------------------------------------------------------*/
void MorphMesh::execute()
{
	std::cout<<"Marking cells"<<std::endl;
	markLockedCells();

	for(auto n : m_mesh->nodes()){
		Node node = m_mesh->get<Node>(n);
		if(node.get<Region>().size() == 4)
			m_surfNodes.push_back(n);
	}

	std::vector<double> durations1;
	std::vector<double> durations2;

	std::cout<<"Start deformation"<<std::endl;
	time_t startglob, endglob;
	time(&startglob);
	for(const auto& p : m_targets){
		std::cout<<"Target point "<<p<<std::endl;
		TCellID n0;
		double min_dist = MAXFLOAT;
		//On cherche le noeud sur le bord du maillage le plus proche, ça sera notre noeud initial
		for(auto n : m_surfNodes){
			math::Point n_point = m_mesh->get<Node>(n).point();
			double dist = n_point.distance(p);
			if(dist < min_dist){
				min_dist = dist;
				n0 = n;
			}
		}
		Node node0 = m_mesh->get<Node>(n0);
		math::Point point0 = node0.point();
		// on calcul le ration d'homothétie pour ce point
		//double H = findHomothetyRatio(p,node0);
		double Y0 = point0.Y();
		double X0 = point0.X();
		double Z = (Y0 - p.Y())/2;
		math::Point tmp;
		computeLocalOrigin(point0, tmp);
		double omega_height0 = Y0 - tmp.Y();
		math::Vector3d t;
		t.setXYZ(0,p.Y()-Y0,0);
		math::Vector3d tz;
		tz.setXYZ(0,0,Z);


		for(auto n : m_mesh->nodes()) {
			if (!m_mesh->isMarked<Node>(n, m_locked)) {
				Node node = m_mesh->get<Node>(n);

				double dist = node.point().distance(point0);
				if (dist <= m_radius) {
					time_t start,end;
					time(&start);
					// Si le point est couvert par la zone à modifier alors on continue le traitement
					//double w = 1 - (dist / m_radius);
					double sigmoid = 1/(1 + exp(-5*( ( (1 - (dist / m_radius) ) *2)-1) ) );

					math::Point n_origin;

					bool origin_locked = computeLocalOrigin(node.point(), n_origin);

					double omega_height = Y0 - n_origin.Y();
					double n_height = node.point().Y() - n_origin.Y();

					math::Vector3d final_t = (n_height / omega_height) * (sigmoid * t);
					math::Point result = node.point() + final_t;

					node.setPoint(result);
					time(&end);
					durations1.emplace_back(double(end - start));
				}

				double distX = fabs(node.X() - X0);
				if (distX <= m_radius) {
					time_t start,end;
					time(&start);
					// double w = 0.5+2*cbrt((1 - (distX / m_radius)) - 0.5);
					double sigmoid = 1 / (1 + exp(-5 * (((1 - (distX / m_radius)) * 2) - 1)));

					math::Point n_origin;
					bool origin_locked = computeLocalOriginZ(node.point(), n_origin);

					double omega_height = fabs(Y0) - fabs(n_origin.Z());
					double n_height = fabs(node.point().Z()) - fabs(n_origin.Z());

					int sens = node.Z() < 0 ? -1 : 1;
					math::Vector3d final_t = (n_height / omega_height) * (sigmoid * (sens * tz));

					math::Point result = node.point() + final_t;
					node.setPoint(result);
					time(&end);
					durations1.emplace_back(double(end - start));
				}
			}
		}
		/*for(auto n : m_mesh->nodes()) {
			if (!m_mesh->isMarked<Node>(n, m_locked)) {
				Node node = m_mesh->get<Node>(n);
				double distX = fabs(node.X() - X0);
				if (distX <= m_radius) {
					// double w = 0.5+2*cbrt((1 - (distX / m_radius)) - 0.5);
					double sigmoid = 1 / (1 + exp(-5 * (((1 - (distX / m_radius)) * 2) - 1)));

					math::Point n_origin;
					bool origin_locked = computeLocalOriginZ(node.point(), n_origin);

					double omega_height = fabs(Y0) - fabs(n_origin.Z());
					double n_height = fabs(node.point().Z()) - fabs(n_origin.Z());

					int sens = node.Z() < 0 ? -1 : 1;
					math::Vector3d final_t = (n_height / omega_height) * (sigmoid * (sens * tz));

					math::Point result = node.point() + final_t;
					node.setPoint(result);
				}
			}
		}*/
	}

	time(&endglob);


	double avg = durations1[0];
	for(int i = 1; i<durations1.size(); i++){
		avg += durations1[i];
	}
	avg /= (double)durations1.size();
	std::cout<<"avg duration for computing local origin "<<avg<<std::setprecision(5);
	std::cout<<" sec for "<<durations1.size()<<" nodes"<<std::endl;
	std::cout<<"globale duration for deformation "<<double(endglob-startglob)<<std::setprecision(5)<<std::endl;

	m_mesh->unmarkAll<Node>(m_locked);
	m_mesh->freeMark<Node>(m_locked);
}
/*----------------------------------------------------------------------------*/
void MorphMesh::markLockedCells()
{
	Variable<int>* locked_regions = m_mesh->newVariable<int,GMDS_REGION>("m_locked");
	std::set<Region> regions;
	std::set<TCellID> nodes;
	nodes.insert(27006);
	m_mesh->mark<Node>(27006,m_locked);
	for (int i = 0; i < 15; ++i) {
		for(auto n : nodes){
			Node node = m_mesh->get<Node>(n);
			for(const auto& r : node.get<Region>()){
				regions.emplace(r);
			}
		}
		nodes.clear();
		for(const auto& r : regions){
			locked_regions->set(r.id(),1);
			for(auto n : r.getIDs<Node>()){
				if(!m_mesh->isMarked<Node>(n,m_locked)){
					nodes.emplace(n);
					m_mesh->mark<Node>(n,m_locked);
				}
			}
		}
		regions.clear();
	}

	for(auto n : m_mesh->nodes()){
		Node node = m_mesh->get<Node>(n);
		if(node.Y() <= 0 ){
			m_mesh->mark<Node>(n,m_locked);
		}
	}

	//Test origine locale décallée
	/*for(auto r : m_mesh->regions()){
		//if((4440<=r && r<=4639) || (5140<=r && r<=5339) || (5840<=r && r<=6039) || (6540<=r && r<=6739)){
			if((37000<=r && r<=37999) || (39000<=r && r<=39999) || (41000<=r && r<=41999) || (43000<=r && r<=43999)){
			locked_regions->set(r,1);
			Region region = m_mesh->get<Region>(r);
			for(auto n : region.getIDs<Node>()){
				m_mesh->mark<Node>(n,m_locked);
			}
		}
	}*/


	for(auto f : m_mesh->faces()){
		Face face = m_mesh->get<Face>(f);

		std::vector<TCellID> f_regions = face.getIDs<Region>();
		if(f_regions.size() == 2){
			//Une face entre un hexa lock et un non lock
			if((locked_regions->value(f_regions[0]) == 1) && (locked_regions->value(f_regions[1]) == 0) ||
			    (locked_regions->value(f_regions[0]) == 0) && (locked_regions->value(f_regions[1]) == 1)){
				locked_faces->set(f,1);
				m_locked_faces.push_back(f);
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
bool MorphMesh::computeLocalOrigin(const math::Point& AP, math::Point& AResult)
{
	/*
	 *  1) on crée un segment entre de AP et (X,Y,0)
	 *  2) On intersecte le segment avec les faces lock
	 *  			- si on a un point alors on retourne ce point
	 *  			- sinon on retourne (X,Y,0)
	 */
	std::vector<math::Point> points;
	points.emplace_back(AP.X(),0,AP.Z());
	points.emplace_back(AP.X()+0.001,0,AP.Z()+0.001);
	points.emplace_back(AP.X()-0.001,0,AP.Z()-0.001);

	for (auto f : m_locked_faces) {
		for(const auto& p : points) {
			//if (locked_faces->value(f) == 1) {
				math::Segment depth(AP, p);
				Face face = m_mesh->get<Face>(f);
				math::Plane plane(face.center(), face.normal());
				// on cherche si la profondeur a une intersection avec un plan d'une face locked
				if (depth.intersect(plane, AResult)) {
					//  si oui on cherche si le point d'intersection est dans la face
					//  On fait ça en créant deux triangles qui couvrent entièrement la face
					//   et qui permettent d'utiliser la méthode isIn()
					std::vector<Node> nodes = face.get<Node>();
					math::Triangle tri1(nodes[0].point(), nodes[1].point(), nodes[2].point());
					math::Triangle tri2(nodes[2].point(), nodes[3].point(), nodes[0].point());
					if (tri1.intersect(depth) || tri2.intersect(depth)) {
						return true;
					}
				}
			//}
		}
	}
	AResult = {AP.X(),0,AP.Z()};

	return false;
}
/*----------------------------------------------------------------------------*/
bool MorphMesh::computeLocalOriginZ(const math::Point& AP, math::Point& AResult)
{
	/*
	 *  1) on crée un segment entre de AP et (X,Y,0)
	 *  2) On intersecte le segment avec les faces lock
	 *  			- si on a un point alors on retourne ce point
	 *  			- sinon on retourne (X,Y,0)
	 */
	std::vector<math::Point> points;
	points.emplace_back(AP.X(),AP.Y(),0);
	points.emplace_back(AP.X()+0.001,AP.Y()+0.001,0);
	points.emplace_back(AP.X()-0.001,AP.Y()-0.001,0);

	for (auto f : m_locked_faces) {
		for(const auto& p : points) {
			//if (locked_faces->value(f) == 1) {
				math::Segment depth(AP,math::Point(AP.X(),AP.Y(),0));
				Face face = m_mesh->get<Face>(f);
				math::Plane plane(face.center(), face.normal());
				// on cherche si la profondeur a une intersection avec un plan d'une face locked
				if (depth.intersect(plane, AResult)) {
					// std::cout<<"test"<<std::endl;
					//  si oui on cherche si le point d'intersection est dans la face
					//  On fait ça en créant deux triangles qui couvrent entièrement la face
					//   et qui permettent d'utiliser la méthode isIn()
					std::vector<Node> nodes = face.get<Node>();
					math::Triangle tri1(nodes[0].point(), nodes[1].point(), nodes[2].point());
					math::Triangle tri2(nodes[2].point(), nodes[3].point(), nodes[0].point());
					if (tri1.intersect(depth) || tri2.intersect(depth)) {
						return true;
					}
				}
			//}
		}
	}

	AResult = {AP.X(),AP.Y(),0};

	return false;
}
/*----------------------------------------------------------------------------*/
double MorphMesh::findHomothetyRatio(const math::Point &AP, const Node& ANode)
{

	math::Point p_node = ANode.point();
	// Le vecteur de déformation initial : La translation du noeud du maillage vers le point
	math::Vector3d deform(AP - ANode.point());


	// On cherche l'origine du noeud initial
	math::Point origin;
	bool origin_locked = computeLocalOrigin(ANode.point(), origin);
	// On fait un vecteur origin->ni pour calculer l'homothétie des deux vecteurs.
	// Faire le vecteur dans ce sens là permet de facilement déduire le ration d'homothétie
	math::Vector3d height(ANode.point() - origin);


	double dist1 = origin.distance(ANode.point());
	double dist2 = origin.distance(AP);

	std::cout<<"dist origin -> node "<<dist1<<std::endl;
	std::cout<<"dist origin -> point "<<dist2<<std::endl;

	return dist2/dist1;
}
/*----------------------------------------------------------------------------*/
}  // namespace MorphMesh
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/