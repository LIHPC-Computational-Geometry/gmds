/*------------------------------------------------------------------------*/
//
// Created by Claire Roche on 21/10/2021.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/Smooth2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
Smooth2D::Smooth2D(Mesh *AMesh,
                   const Variable<int>* AVarBnd,
                   const int ANbIterations) {
	m_mesh = AMesh;
	m_node_constrained = AVarBnd;
	m_nb_max_iterations = ANbIterations;
}
/*------------------------------------------------------------------------*/
void Smooth2D::setNbIterations(const int ANbIterations)
{
	m_nb_max_iterations=ANbIterations;
}
/*------------------------------------------------------------------------*/
Smooth2D::STATUS Smooth2D::execute()
{
    //we traverse all the nodes of the whole mesh to get the nodes that
    //are not constrained, i.e are so free
	for(auto n_id:m_mesh->nodes()){
		if(m_node_constrained->value(n_id)!=1){
			m_free_nodes.push_back(n_id);
		}
	}
    //for all free nodes, we build and store stencils locally
    buildStencils();


	 /*
	 std::cout << "TEST 1 :" << std::endl ;
	 math::Point A(0,0,0);
	 math::Point B(1,0,0);
	 math::Point C(3,0,0);
	 math::Point Test_Mid ;
	 Test_Mid = FindMidBranche(A, B, C) ;

	 std::cout << "----------------------" << std::endl ;
	 std::cout << "TEST 2 :" << std::endl ;
	 A = (0,0,0);
	 B = (0,1,0);
	 B.X() = 0,
	 B.Y() = 1;
	 B.Z() = 0;
	 C.X() = 2;
	 C.Y() = 1;
	 C.Z() = 0;
	 Test_Mid = FindMidBranche(A, B, C) ;
	  */

	    return Smooth2D::SUCCESS;
}
/*------------------------------------------------------------------------*/
void Smooth2D::buildStencils() {
    for(auto n_id:m_free_nodes){
        Node current_node = m_mesh->get<Node>(n_id);
        std::vector<Face> current_faces = current_node.get<Face>();
        // FIRST FACE
        Face f0 = current_faces[0];
        Node f0_n0;
        Node f0_n1;
        f0.getAdjacentNodes(current_node,f0_n0,f0_n1);
        stencil current_stencil;
        current_stencil.val[1][1] = current_node.id();
        current_stencil.val[1][0] = f0_n0.id();
        current_stencil.val[0][1] = f0_n1.id();
        std::vector<TCellID> f_nodes=f0.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f0_n0.id() &&
               f_n_id!=f0_n1.id()) {
                current_stencil.val[0][0] = f_n_id;
            }
        }
        //SECOND FACE
        std::vector<TCellID> faces_11_10 = m_mesh->getCommonFaces(current_node,f0_n0);
        TCellID f1_id = NullID;
        if(faces_11_10[0]==f0.id()){
            f1_id=faces_11_10[1];
        }
        else{
            f1_id=faces_11_10[0];
        }
        Face f1 = m_mesh->get<Face>(f1_id);
        Node f1_n0;
        Node f1_n1;
        f1.getAdjacentNodes(current_node,f1_n0,f1_n1);
        if(f1_n0.id()==current_stencil.val[1][0])
            current_stencil.val[2][1]=f1_n1.id();
        else
            current_stencil.val[2][1]=f1_n0.id();

        f_nodes=f1.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f1_n0.id() &&
               f_n_id!=f1_n1.id()) {
                current_stencil.val[2][0] = f_n_id;
            }
        }
        //THIRD FACE
        std::vector<TCellID> faces_11_21 =
                m_mesh->getCommonFaces(current_node,
                                       m_mesh->get<Node>(current_stencil.val[2][1]));
        TCellID f2_id = NullID;
        if(faces_11_21[0]==f1.id()){
            f2_id=faces_11_21[1];
        }
        else{
            f2_id=faces_11_21[0];
        }
        Face f2 = m_mesh->get<Face>(f2_id);
        Node f2_n0;
        Node f2_n1;
        f2.getAdjacentNodes(current_node,f2_n0,f2_n1);

        if(f2_n0.id()==current_stencil.val[2][1])
            current_stencil.val[1][2]=f2_n1.id();
        else
            current_stencil.val[1][2]=f2_n0.id();

        f_nodes=f2.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f2_n0.id() &&
               f_n_id!=f2_n1.id()) {
                current_stencil.val[2][2] = f_n_id;
            }
        }

        //FOURTH FACE
        std::vector<TCellID> faces_11_12 =
                m_mesh->getCommonFaces(current_node,
                                       m_mesh->get<Node>(current_stencil.val[1][2]));
        TCellID f3_id = NullID;
        if(faces_11_12[0]==f2.id()){
            f3_id=faces_11_12[1];
        }
        else{
            f3_id=faces_11_12[0];
        }
        Face f3 = m_mesh->get<Face>(f3_id);

        f_nodes=f3.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=current_stencil.val[1][2] &&
               f_n_id!=current_stencil.val[0][1]) {
                current_stencil.val[0][2] = f_n_id;
            }
        }

        m_stencil.insert(std::make_pair(n_id, current_stencil));
		  // Affichage du stencil pour chaque noeud intérieur
		  /*
        std::cout<<"Node "<<n_id<<std::endl;
        for(auto i=0;i<3;i++){
            for(auto j=0;j<3;j++){
                std::cout << m_stencil[n_id].val[i][j] << " ";
            }
            std::cout<<std::endl;
        }
        */
    }
}
/*------------------------------------------------------------------------*/


/*
void Smooth2D::PerturbationMaillage(const Variable<int>* var_bnd, const double dx, const double dy) {
	// ------------ PERTURBATION ALEATOIRE DU MAILLAGE ------------
	// Initialisation de la variable n qui est un noeud et n_coords qui récupère les coordonnées du noeud
	Node n = m_mesh->get<Node>(1);
	math::Point n_coords = n.point();
	double c1 = 0.0;
	double c2 = 0.0;
	for(auto id:m_mesh->nodes()){
		if(var_bnd->value(id)==0){
			n = m_mesh->get<Node>(id);
			n_coords = n.point();
			c1 = rand()/ double(RAND_MAX) ;
			c2 = rand()/ double(RAND_MAX) ;
			//std::cout << "Numéro node :" << n << std::endl ;
			//std::cout << "Nombre aléatoire c1 =" << c1 << std::endl ;
			n.setX( n_coords.X() + (2.0*c1-1.0) * dx);
			n.setY( n_coords.Y() + (2.0*c2-1.0) * dy);
		}
	}
	// ---------------------------------------------------------------
}
 */





/*------------------------------------------------------------------------*/
// Fonction FindMidBranche : Si on considère une branche composée de 3 points, cette fonction retourne le point
//                            positionné au milieu de cette branche
// En entrée : A, B, C -> 3 points
// En sortie : le point milieu
math::Point Smooth2D::FindMidBranche(const math::Point A, const math::Point B, const math::Point C) {
	math::Point Point_Milieu ;
	double norme_1 = sqrt( pow((A.X()-B.X()),2) + pow(A.Y()-B.Y(),2) ) ;
	double norme_2 = sqrt( pow((B.X()-C.X()),2) + pow(B.Y()-C.Y(),2) ) ;
	double norme_branche = norme_1 + norme_2 ;
	double norme_milieu = norme_branche / 2.0 ;

	if (norme_milieu <= norme_1){
		math::Vector3d Vec( B.X()-A.X(), B.Y()-A.Y(), 0 ) ;
		Vec.normalize();
		Point_Milieu = A + norme_milieu*Vec ;
	}
	else if (norme_milieu > norme_1){
		math::Vector3d Vec( B.X()-C.X(), B.Y()-C.Y(), 0 ) ;
		Vec.normalize();
		Point_Milieu = C + norme_milieu*Vec ;
	}

	std::cout << "Point milieu :" << Point_Milieu << std::endl ;
	return Point_Milieu;
}
/*------------------------------------------------------------------------*/