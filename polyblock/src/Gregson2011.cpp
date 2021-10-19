/*---------------------------------------------------------------------------*/
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <ctime>
/*---------------------------------------------------------------------------*/
#include <gmds/math/Vector.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/VTKWriter.h>
/*---------------------------------------------------------------------------*/
#include <gmds/polyblock/Gregson2011.h>
#include <gmds/io/IGMeshIOService.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace Eigen;
/*---------------------------------------------------------------------------*/
typedef Matrix<double,3,3> Matrix3;
typedef Matrix <double, 2, 1> vect2D;
typedef SparseMatrix<double> SpMat;
/*---------------------------------------------------------------------------*/
//																GREGSON2011_2D
/*---------------------------------------------------------------------------*/
Gregson2011_2D::Gregson2011_2D(Mesh *mesh){
    this->mesh = mesh;
    markOnBoundary();
}
/*---------------------------------------------------------------------------*/
Mesh Gregson2011_2D::getPolycube(){
    Mesh returnMesh(this->mesh->getModel());
    return returnMesh;
}
/*---------------------------------------------------------------------------*/
math::Vector eigen2math(vect2D v){
    return math::Vector(v.x(),v.y(),0);
}
/*---------------------------------------------------------------------------*/
Mesh Gregson2011_2D::getMesh(){
    return *(this->mesh);
}
/*---------------------------------------------------------------------------*/
void Gregson2011_2D::markOnBoundary(){
	mark_smooth_bnd_nodes = mesh->newMark<Node>();
	mark_bnd_nodes = mesh->newMark<Node>();
	mark_bnd_egdes = mesh->newMark<Edge>();

    for(auto e_id: mesh->edges())
    {
        Edge e = mesh->get<Edge>(e_id);
        if (e.nbFaces() == 1) {
            mesh->mark(e,mark_bnd_egdes);
            std::vector<TCellID> nodes = e.getIDs<Node>();
            for (auto n:nodes){
                mesh->mark(mesh->get<Node>(n),mark_bnd_nodes);
            }
        }
    }

    // a node is smooth if the dot product between the two normal of its boudary edge are > 0.5
    for (auto n_id : mesh->nodes())
    {
        Node n = mesh->get<Node>(n_id);
        if (mesh->isMarked(n,mark_bnd_nodes)){
            std::vector<TCellID> edges = n.getIDs<Edge>();
            std::vector<TCellID> bnd_edge;
            for (auto e:edges){
                if (mesh->isMarked(mesh->get<Edge>(e),mark_bnd_egdes)){
                    bnd_edge.push_back(e);
                }
            }

            if (normal(mesh->get<Edge>(bnd_edge[0])).transpose() * normal(mesh->get<Edge>(bnd_edge[1])) > 0.8){
                mesh->mark(n,mark_smooth_bnd_nodes);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void Gregson2011_2D::updateBoundary(){
	for (auto n_id: mesh->nodes()) {
		Node n = mesh->get<Node>(n_id);
		if (mesh->isMarked(n,mark_bnd_nodes)){
			std::vector<TCellID> edges = n.getIDs<Edge>();
			std::vector<TCellID> bnd_edge;
			for (auto e:edges){
				if (mesh->isMarked(mesh->get<Edge>(e),mark_bnd_egdes)){
					bnd_edge.push_back(e);
				}
			}
			if (normal(mesh->get<Edge>(bnd_edge[0])).transpose() * normal(mesh->get<Edge>(bnd_edge[1])) > 0.8){
				mesh->mark(n,mark_smooth_bnd_nodes);
			} else if (mesh->isMarked(n,mark_smooth_bnd_nodes)) {
				mesh->unmark(n,mark_smooth_bnd_nodes);
			}

		}
	}
}

/*----------------------------------------------------------------------------*/
void Gregson2011_2D::transformMeshToPolycube(){
    
	// // wanted variable
	Variable<int>* edge_chosenNormal_ID = mesh->newVariable<int,GMDS_EDGE>("ChosenNormal_ID");
	Variable<vect2D>* edge_chosenNormal = mesh->newVariable<vect2D,GMDS_EDGE>("ChosenNormal");

	//display variable
	Variable<math::Vector>* display_normals = mesh->newVariable<math::Vector,GMDS_NODE>("d_normal");
	Variable<math::Vector>* display_chosen_normals = mesh->newVariable<math::Vector,GMDS_NODE>("d_ChosenNormal");
	Variable<math::Vector>* rotated_normal = mesh->newVariable<math::Vector,GMDS_NODE>("rotated_normal");

	//usefull running information
	Variable<vect2D>* point_normals = mesh->newVariable<vect2D,GMDS_NODE>("normal");
	Variable<double>* point_angle2chosen = mesh->newVariable<double,GMDS_NODE>("point_angle2chosen");
	Variable<vect2D>* edge_normals = mesh->newVariable<vect2D,GMDS_EDGE>("normal");



	Matrix<double,Dynamic,Dynamic> adjacency = Matrix<double,Dynamic,Dynamic>::Zero(mesh->getNbNodes(),mesh->getNbNodes());
	// calcul la matrice d'adjacence
	for(auto e_id : mesh->edges()){
		Edge e = mesh->get<Edge>(e_id);
		std::vector<TCellID> nodes = e.getIDs<Node>();
		adjacency(nodes[0],nodes[1]) = 1;
		adjacency(nodes[1],nodes[0]) = 1;
	}
	for (auto n_id : mesh->nodes()) {
		Node n = mesh->get<Node>(n_id);
		adjacency(n.id(),n.id()) = -n.nbEdges();

	}
    Matrix<double,Dynamic,Dynamic> poisson = - adjacency;

	Matrix<double,Dynamic,1> smooth_bnd_angle(mesh->getNbNodes(),1);
	smooth_bnd_angle<<  Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
	Matrix<double,Dynamic,1> generalised_angle(mesh->getNbNodes(),1);
	Matrix<double,Dynamic,2> new_position(mesh->getNbNodes(),2);
	Matrix<double,Dynamic,2> poisson_2nd_member(mesh->getNbNodes(),2);


	Variable<int>* on_bnd = mesh->newVariable<int,GMDS_NODE>("on_bnd");



	for (int iteration = 0; iteration < 3; iteration++) {
		std::cout << "---------" << "Iteration " << iteration+1 << "---------" << '\n';
		for(auto e_id : mesh->edges())
		{
			Edge e= mesh->get<Edge>(e_id);
			(*edge_normals)[e.id()] = normal(e);
			(*edge_chosenNormal)[e.id()] = getClosestCanonical((*edge_normals)[e.id()]);
			(*edge_chosenNormal_ID)[e.id()] = getCanonicalId((*edge_chosenNormal)[e.id()]);
		}


		//calcul des vecteurs normal aux points
        for (auto n_id : mesh->nodes()) {
            Node n = mesh->get<Node>(n_id);
            if (mesh->isMarked(n,mark_smooth_bnd_nodes)){
				std::vector<TCellID> edges = n.getIDs<Edge>();
				(*point_normals)[n.id()]=vect2D(0.,0.);
                for (auto e:edges){
                    if (mesh->isMarked(mesh->get<Edge>(e), mark_bnd_egdes)){
                        (*point_normals)[n.id()] += (*edge_normals)[e];
                    }
                }
                (*point_normals)[n.id()].normalize();

                (*display_normals)[n.id()] = eigen2math((*point_normals)[n.id()]);
                (*display_chosen_normals)[n.id()] = eigen2math(getClosestCanonical((*point_normals)[n.id()]) );
                vect2D v = (*point_normals)[n.id()];
                (*point_angle2chosen)[n.id()] = ((vect2D(v.y(),-v.x()).transpose() * getClosestCanonical((*point_normals)[n.id()]) > 0 ) ? -1 : 1)
                                                *acos(v.transpose() * getClosestCanonical((*point_normals)[n.id()]));
                smooth_bnd_angle(n.id(),0)  =  (*point_angle2chosen)[n.id()];
                adjacency.row(n.id()) = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
                (*on_bnd)[n.id()] = 1;
                adjacency(n.id(),n.id()) =1;
            }
        }

        // adjacency.row(n.id()) = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
		// adjacency(n.id(),n.id()) =1;
        // std::cout << adjacency << '\n';
		// std::cout << "smooth" << endl <<smooth_bnd_angle << '\n';
		generalised_angle = adjacency.fullPivLu().solve(smooth_bnd_angle);
		// std::cout << "general" << endl <<generalised_angle << '\n';

		double angle1,angle2;
        for (auto n_id : mesh->nodes()) {
            Node n = mesh->get<Node>(n_id);
            poisson_2nd_member.row(n.id()) << 0 ,0;
			std::vector<TCellID> edges = n.getIDs<Edge>();
			angle1 = generalised_angle(n.id());
			for (auto e:edges){
				std::vector<TCellID> nodes = (mesh->get<Edge>(e)).getIDs<Node>();
				for (auto e_n : nodes) {
                    if (!(e_n == n.id())){
                        angle2 = generalised_angle(e_n);
                        Matrix2d m;
                        m << cos(angle1)+cos(angle2), -sin(angle1)-sin(angle2),
                                sin(angle1)+sin(angle2), cos(angle1)+cos(angle2);
                        // std::cout << m << '\n';
                        poisson_2nd_member.row(n.id()) += m * vect2D(n.X() - (mesh->get<Node>(e_n)).X(),
                                                                     n.Y() - (mesh->get<Node>(e_n)).Y())/2;

                    }
                }
			}
		}


		poisson.row(0) = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
		poisson(0,0) =1;
		poisson_2nd_member.row(0) = vect2D((mesh->get<Node>(0)).X(), (mesh->get<Node>(0)).Y());

		// std::cout << "Poisson = " << endl << poisson  << '\n';
		// std::cout << "poisson_2nd_member = " << endl << poisson_2nd_member  << '\n';

		new_position = poisson.fullPivLu().solve(poisson_2nd_member);
		//
		// std::cout << "new_position = "<<endl << new_position << '\n';
		// std::cout << "result = " << endl << poisson * new_position <<'\n';


        for (auto n_id : mesh->nodes()) {
            Node n = mesh->get<Node>(n_id);
            // std::cout << "en X  :: "<<n.X() <<" devient "  << new_position(n.id(),0) << '\n';
			// std::cout << "en Y  :: "<<n.Y() <<" devient "  << new_position(n.id(),1) << '\n';
			n.setX(new_position(n.id(),0));
			n.setY(new_position(n.id(),1));
		}

        
		//display the rotation planned, without moving the point
        for (auto n_id : mesh->nodes()) {
            Node n = mesh->get<Node>(n_id);
            if (mesh->isMarked(n,mark_smooth_bnd_nodes)){
				std::vector<TCellID> edges = n.getIDs<Edge>();
				(*rotated_normal)[n.id()]=math::Vector(0.,0.,0.);
				for (auto e:edges){
					if (mesh->isMarked(mesh->get<Edge>(e), mark_bnd_egdes)){
						std::vector<TCellID> nodes = (mesh->get<Edge>(e)).getIDs<Node>();
						vect2D rotated_edge_normal;
						double cond = 0;
						if (nodes[0]!= n.id()){
							rotated_edge_normal = vect2D(new_position(nodes[0],1) - new_position(n.id(),1),new_position(n.id(),0) - new_position(nodes[0],0));
							rotated_edge_normal.normalize();
							cond =  ((rotated_edge_normal.transpose() * getClosestCanonical((*edge_normals)[e]) < 0 ) ? -1 : 1);
							rotated_edge_normal *= cond;
						}
						else {
							rotated_edge_normal = vect2D(new_position(nodes[1],1) - new_position(n.id(),1),new_position(n.id(),0) - new_position(nodes[1],0));
							rotated_edge_normal.normalize();
							cond =  ((rotated_edge_normal.transpose() * getClosestCanonical((*edge_normals)[e]) < 0 ) ? -1 : 1);
							rotated_edge_normal *= cond;
						}
						(*rotated_normal)[n.id()] =(*rotated_normal)[n.id()] + eigen2math(rotated_edge_normal);
					}
				}
				(*rotated_normal)[n.id()].normalize();
			}

		}
		this->updateBoundary();

        gmds::IGMeshIOService ioService(mesh);
        gmds::VTKWriter w_2D(&ioService);
        w_2D.setCellOptions(gmds::N|gmds::E|gmds::F);
        w_2D.setDataOptions(gmds::N|gmds::E|gmds::F);
 		w_2D.write("result2D_"+to_string(iteration+1));

		// modifier la position des points avec le system de poisson

		// rebouclé

	}
    for (auto e_id : mesh->edges()) {
        Edge e = mesh->get<Edge>(e_id);
        (*edge_normals)         [e_id] = normal(e);
		(*edge_chosenNormal)    [e_id] = getClosestCanonical((*edge_normals)[e_id]);
		(*edge_chosenNormal_ID) [e_id] = getCanonicalId((*edge_chosenNormal)[e_id]);
	}
    for (auto n_id : mesh->nodes()) {
        Node n = mesh->get<Node>(n_id);
		if (mesh->isMarked(n,mark_smooth_bnd_nodes)){
            std::vector<TCellID> edges = n.getIDs<Edge>();
            (*point_normals)[n.id()]=vect2D(0.,0.);
            for (auto e:edges){
                if (mesh->isMarked(mesh->get<Edge>(e), mark_bnd_egdes)){
                    (*point_normals)[n.id()] += (*edge_normals)[e];
                }
            }
            (*point_normals)[n.id()].normalize();

            (*display_normals)[n.id()] = eigen2math((*point_normals)[n.id()]);
            (*display_chosen_normals)[n.id()] = eigen2math(getClosestCanonical((*point_normals)[n.id()]) );
            (*point_angle2chosen)[n.id()] = acos((*point_normals)[n.id()].transpose()
                                                 * getClosestCanonical((*point_normals)[n.id()]));
			smooth_bnd_angle(n.id(),0) =(*point_angle2chosen)[n.id()];
			adjacency.row(n.id()) = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
			adjacency(n.id(),n.id()) =1;
		}
	}

    this->projectOnPolycube(edge_chosenNormal);

}
/*---------------------------------------------------------------------------*/
void Gregson2011_2D::projectOnPolycube(Variable<vect2D>* edge_chosenNormal){
	// We are now given a closest mesh to a polycube, we divide boundary into chart
	typedef std::vector<TCellID> chart;
	std::vector<chart> charts;
	int current_chart = 0;
	chart edge_to_explore;
	int marked_edge = mesh->newMark<Edge>();
	for (auto e_id: mesh->edges()) {
        Edge source_e= mesh->get<Edge>(e_id);
        if (!(mesh->isMarked(source_e,marked_edge)) && (mesh->isMarked(source_e,mark_bnd_egdes))){
            charts.push_back(chart());
            edge_to_explore.push_back(source_e.id());
            while (!edge_to_explore.empty()) {
                source_e = mesh->get<Edge>(edge_to_explore.back());
                edge_to_explore.pop_back();
                charts[current_chart].push_back(source_e.id());
                mesh->mark(source_e,marked_edge);
                std::vector<TCellID> nodes = source_e.getIDs<Node>();
                std::vector<TCellID> edges = (mesh->get<Node>(nodes[0])).getIDs<Edge>();
                for (auto e:edges){
                    if ( (*edge_chosenNormal)[source_e.id()] == (*edge_chosenNormal)[e]){
                        if (!mesh->isMarked(mesh->get<Edge>(e),marked_edge)&& mesh->isMarked(mesh->get<Edge>(e),mark_bnd_egdes)){
                            edge_to_explore.push_back(e);
                        }
                    }
                }
                edges = (mesh->get<Node>(nodes[1])).getIDs<Edge>();
                for (auto e:edges){
                    if ((*edge_chosenNormal)[source_e.id()] == (*edge_chosenNormal)[e]){
                        if (!mesh->isMarked(mesh->get<Edge>(e),marked_edge)&& mesh->isMarked(mesh->get<Edge>(e),mark_bnd_egdes)){
							edge_to_explore.push_back(e);
						}
					}
				}
			}
			current_chart++;
		}
	}

	mesh->unmarkAll<Edge>(marked_edge);
	mesh->freeMark<Edge>(marked_edge);

	// pushing into polycube, better commented for testing other part
	Variable<int>* node_updatedX = mesh->newVariable<int,GMDS_NODE>("updatedX");
	Variable<int>* node_updatedY = mesh->newVariable<int,GMDS_NODE>("updatedY");

	Matrix<double,Dynamic,Dynamic> adjacency = Matrix<double,Dynamic,Dynamic>::Zero(mesh->getNbNodes(),mesh->getNbNodes());

	for(auto e_id : mesh->edges()){
		Edge e = mesh->get<Edge>(e_id);
		std::vector<TCellID> nodes = e.getIDs<Node>();
		adjacency(nodes[0],nodes[1]) = -1;
		adjacency(nodes[1],nodes[0]) = -1;
	}
	for (auto n_id:mesh->nodes()) {
		Node n = mesh->get<Node>(n_id);
		adjacency(n.id(),n.id()) = n.nbEdges();
	}
	Matrix<double,Dynamic,2> scnd_member = Matrix<double,Dynamic,2>::Zero(mesh->getNbNodes(),2);
	Matrix<double,Dynamic,2> new_position ;

    for (auto n_id:mesh->nodes()) {
        Node n = mesh->get<Node>(n_id);
        if (mesh->isMarked(n,mark_bnd_nodes)){
			scnd_member.row(n.id()) = vect2D(n.X(),n.Y());
		}
	}

	vect2D working_vector;
	std::vector<TCellID> chart_nodes;
	double max_value,set_value;
	// Now we want to push each chart to a border of a polycube
	for (auto c : charts){
		max_value =-10E80;
		set_value =0;
		for (auto e :c){
			chart_nodes = (mesh->get<Edge>(e)).getIDs<Node>();
			for (auto node:chart_nodes){
				Node n = mesh->get<Node>(node);
				working_vector = vect2D(n.X(),n.Y());
				if (working_vector.transpose()*	(*edge_chosenNormal)[e] > max_value){
					set_value = ((vect2D(1.,1.).transpose() *(*edge_chosenNormal)[e] )*(working_vector.transpose())*	(*edge_chosenNormal)[e] );
					max_value = working_vector.transpose()*	(*edge_chosenNormal)[e];
				}
			}
		}
		for (auto e :c){
			chart_nodes = (mesh->get<Edge>(e)).getIDs<Node>();
			if ((*edge_chosenNormal)[e].x() != 0 ){
				for (auto node:chart_nodes){
					// (mesh->get<Node>(node)).setX(set_value); // a decommenter pour comparer avec projection non etalée
					scnd_member(node,0)=set_value;
					adjacency.row(node) = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
					adjacency(node,node)=1;
					(*node_updatedX)[node] = 1;
				}
			}
			else
			{
				for (auto node:chart_nodes){
					// (mesh->get<Node>(node)).setY(set_value);
					scnd_member(node,1)=set_value;
					adjacency.row(node) = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
					adjacency(node,node)=1;
					(*node_updatedY)[node] = 1;
				}
			}
		}
	}

	new_position = adjacency.fullPivLu().solve(scnd_member);

    for (auto n_id:mesh->nodes()) {
        Node n = mesh->get<Node>(n_id);
        n.setX(new_position(n.id(),0));
		n.setY(new_position(n.id(),1));
	}

}
/*----------------------------------------------------------------------------*/
vect2D Gregson2011_2D::normal(Edge e)
{
    // we extract the two nodes of the edge
    std::vector<TCellID> Nodes = e.getIDs<Node>();
    // deduce the vector of the edge
    vect2D edge_vector((mesh->get<Node>(Nodes[0])).X() - (mesh->get<Node>(Nodes[1])).X(),
                       (mesh->get<Node>(Nodes[0])).Y() - (mesh->get<Node>(Nodes[1])).Y() );
    // compute a given normal
    vect2D normal(-edge_vector.y(),edge_vector.x());

    // we now need to see if the normal goes to the opposite side of the mesh (goes out of the mesh)
    // we take the node of the face
    std::vector<TCellID> f_nodes = (mesh->get<Face>(e.getIDs<Face>()[0])).getIDs<Node>();
    vect2D direction_comparator;
    // we locate the node which is not on the face
    for (auto n:f_nodes){
        if ((n != Nodes[0]) && (n != Nodes[1])){
            // we compute the vector from a point of the edge to the extra point
            direction_comparator = vect2D( (mesh->get<Node>(Nodes[0])).X() - (mesh->get<Node>(n)).X(),
                                           (mesh->get<Node>(Nodes[0])).Y() - (mesh->get<Node>(n)).Y() );
        }
    }
    // if the dot product is positive, the two vector go to the general same direction
    if (normal.transpose() * direction_comparator < 0){
        // and we want the normal going the other way around
        normal = -normal;
    }
    normal.normalize();
    return normal;
}
/*---------------------------------------------------------------------------*/
vect2D Gregson2011_2D::getClosestCanonical(vect2D v){ // pas beau
    vect2D closest = vect2D(1.,0.);
	double dotproduct = closest.transpose() * v;
	if (vect2D(-1.,0.).transpose() * v > dotproduct+0.05){
		dotproduct = vect2D(-1.,0.).transpose() * v;
		closest = vect2D(-1.,0.);
	}
	if (vect2D(0.,1.).transpose() * v > dotproduct+0.05){
		dotproduct = vect2D(0.,1.).transpose() * v;
		closest = vect2D(0.,1.);
	}
	if (vect2D(0.,-1.).transpose() * v > dotproduct+0.05){
		dotproduct = vect2D(0.,-1.).transpose() * v;
		closest = vect2D(0.,-1.);
	}
	return closest;
}
/*---------------------------------------------------------------------------*/
int Gregson2011_2D::getCanonicalId(vect2D v){ //pas tres beau
	if (v == vect2D(1.,0.))
		return 0;
		else if (v == vect2D(-1.,0.))
			return 1;
			else if (v == vect2D(0.,1.))
			 	return 2;
				else if (v == vect2D(0.,-1.))
					return 3;
	std::cout << "Non canonical vector passed throught getCanonicalId" << v
	<< "returning 4" <<'\n';
	return 4;
}
/*---------------------------------------------------------------------------*/
