/*----------------------------------------------------------------------------*/
#include "gmds/morphMesh/MorphMesh.h"
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
}
/*----------------------------------------------------------------------------*/
MorphMesh::~MorphMesh() {}
/*----------------------------------------------------------------------------*/
void
MorphMesh::execute()
{
	markLockedCells();

	for(auto n : m_mesh->nodes()){
		Node node = m_mesh->get<Node>(n);
		if(node.get<Region>().size() == 4)
			m_surfNodes.push_back(n);
	}

	for(const auto& p : m_targets){
		std::cout<<"Traitement point "<<p<<std::endl;
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
		double H = findHomothetyRatio(p,node0);


		for(auto n : m_mesh->nodes()){
			if(!m_mesh->isMarked<Node>(n,m_locked)) {
				Node node = m_mesh->get<Node>(n);
				math::Point verifPoint = node.point();
				verifPoint.setY(point0.Y());     // On fait ça pour mettre les deux points sur le même plan

				double dist = verifPoint.distance(point0);
				if (dist <= m_radius) {
					std::cout << "Dans le radius" << std::endl;
					// Si le point est couvert par la zone à modifier alors on continue le traitement
					double w = dist / m_radius;
					std::cout << "Poids = " << w << std::endl;

					math::Point n_origin = computeLocalOrigin(node.point());
					math::Vector3d height = math::Vector3d(node.point() - n_origin);

					std::cout << "before H = " << H << std::endl;
					double localH = H + ((1 - H) * w);     // On réduit le ratio d'homothétie par le poids
					std::cout << "after H = " << localH << std::endl;

					math::Point result = n_origin + height * localH;

					std::cout << "before " << node << std::endl;
					node.setPoint(result);
					std::cout << "after " << node << std::endl;
				}
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void
MorphMesh::markLockedCells()
{
	/** ==================================
	 * 				 !! WIP !!
	    =================================

	    1) this method should also lock nodes given by user
	    2) maybe we should lock cells too, edge, faces or region
	*/
	for(auto n : m_mesh->nodes()){
		Node node = m_mesh->get<Node>(n);
		if(node.Y() <= 0 ){
			m_mesh->mark<Node>(n,m_locked);
		}
	}
}
/*----------------------------------------------------------------------------*/
math::Point
MorphMesh::computeLocalOrigin(const math::Point& AP)
{
	/*
	 *  1) on crée un segment entre de AP et (X,Y,0)
	 *  2) On intersecte le segment avec les faces lock
	 *  			- si on a un point alors on retourne ce point
	 *  			- sinon on retourne (X,Y,0)
	 */

	return {AP.X(),0,AP.Z()};
}
/*----------------------------------------------------------------------------*/
double MorphMesh::findHomothetyRatio(const math::Point &AP, const Node& ANode)
{

	math::Point p_node = ANode.point();
	// Le vecteur de déformation initial : La translation du noeud du maillage vers le point
	math::Vector3d deform(AP - ANode.point());


	// On cherche l'origine du noeud initial
	math::Point origin = computeLocalOrigin(ANode.point());
	// On fait un vecteur origin->ni pour calculer l'homothétie des deux vecteurs.
	// Faire le vecteur dans ce sens là permet de facilement déduire le ration d'homothétie
	math::Vector3d height(ANode.point() - origin);


	double dist1 = origin.distance(ANode.point());
	double dist2 = origin.distance(AP);


	return dist2/dist1;
}
/*----------------------------------------------------------------------------*/
}  // namespace MorphMesh
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/