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

	 math::Point H1, H2, H3;
	 math::Point V1, V2, V3;
	 math::Point A, B, C;
	 Node noeud_voisin;
	 for(auto n_id:m_free_nodes){
		 std::cout << "Noeud :" << n_id << std::endl ;
		 // Noeuds de la première branche (verticale)
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][0]);
		 A = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][0]);
		 B = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][0]);
		 C = noeud_voisin.point();
		 V1 = FindMidBranche(A, B, C);

		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][1]);
		 A = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][1]);
		 B = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][1]);
		 C = noeud_voisin.point();
		 V2 = FindMidBranche(A, B, C);

		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][2]);
		 A = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][2]);
		 B = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][2]);
		 C = noeud_voisin.point();
		 V3 = FindMidBranche(A, B, C);

		 // Noeuds de la seconde branche (horizontale)
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][0]);
		 A = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][1]);
		 B = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][2]);
		 C = noeud_voisin.point();
		 H1 = FindMidBranche(A, B, C);

		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][0]);
		 A = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][1]);
		 B = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][2]);
		 C = noeud_voisin.point();
		 H2 = FindMidBranche(A, B, C);

		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][0]);
		 A = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][1]);
		 B = noeud_voisin.point();
		 noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][2]);
		 C = noeud_voisin.point();
		 H3 = FindMidBranche(A, B, C);
	 }

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

	 bool intersection ;
	 std::cout << " TESTS D'INTERSECTIONS " << std::endl ;
	 std::cout << "TEST 1 :" << std::endl ;
	 math::Point Test_A(0,0,0);
	 math::Point Test_B(2,0,0);
	 math::Point Test_C(1,1,0);
	 math::Point Test_D(1,-1,0);
	 intersection = IntersectionSegments2D( Test_A, Test_B, Test_C, Test_D) ;
	 std::cout << "Intersection ?" << intersection << std::endl ;

	 // Intersection non trouvée (et pourtant, elle existe)
	 std::cout << "TEST 2 :" << std::endl ;
	 Test_A.X() = 0;
	 Test_A.Y() = 0;
	 Test_A.Z() = 0;
	 Test_B.X() = 2,
	 Test_B.Y() = 0;
	 Test_B.Z() = 0;
	 Test_C.X() = 2;
	 Test_C.Y() = 0;
	 Test_C.Z() = 0;
	 Test_D.X() = 2;
	 Test_D.Y() = -2;
	 Test_D.Z() = 0;

	 std::cout << "Point D" << Test_D << std::endl ;
	 intersection = IntersectionSegments2D( Test_A, Test_B, Test_C, Test_D) ;
	 std::cout << "Intersection ?" << intersection << std::endl ;




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

	/*
	std::cout << "Point A :" << A << std::endl;
	std::cout << "Point B :" << B << std::endl;
	std::cout << "Point C :" << C << std::endl;
	 */
	std::cout << "Point milieu :" << Point_Milieu << std::endl ;
	std::cout << "--------------------------------------" << std::endl;
	return Point_Milieu;
}
/*------------------------------------------------------------------------*/







/*------------------------------------------------------------------------*/
// Fonction IntersectionSegments2D : Regarde si deux segments se coupent
// En entrée : A, B, C, D -> 4 points, A et B forment le premier segment, C et D forment le deuxième segment
// En sortie : intersection -> bool, false si il n'y a pas d'intersection, true sinon
bool Smooth2D::IntersectionSegments2D(const math::Point A, const math::Point B, const math::Point C, const math::Point D) {
	bool intersection(false);
	math::Vector3d Vect_1( B.X()-A.X(), B.Y()-A.Y(), 0 );
	math::Vector3d Vect_2( D.X()-C.X(), D.Y()-C.Y(), 0 );
	Vect_1.normalize();
	Vect_2.normalize();
	if ( (Vect_1 != Vect_2) && (Vect_1 != -Vect_2)){
		math::Point Point_intersection(0,0,0);
		bool M_sur_seg1(false);
		bool M_sur_seg2(false);
		// Calcul de l'intersection des deux droites
		// Résolution d'un système 2x2
		// ATTENTION : ERREUR DANS LE CALCUL DU POINT D'INTERSECTION
		Point_intersection.X() = A.X() + ( -Vect_2.Y()*( C.X()-A.X() ) + Vect_2.X()*( C.Y()-A.Y() ) )*Vect_1.X() ;
		Point_intersection.Y() = A.Y() + ( -Vect_2.Y()*( C.X()-A.X() ) + Vect_2.X()*( C.Y()-A.Y() ) )*Vect_1.Y() ;
		std::cout << "HELLO" << std::endl ;
		std::cout << "M :" << Point_intersection << std::endl ;
		// Est-ce que le point d'intersection appartient aux deux segments ?
		// On note que le point d'intersection est appelé M dans la suite pour alléger les notations
		// -> Test segment 1 [A,B]
		double norme_AB = sqrt( pow((A.X()-B.X()),2) + pow(A.Y()-B.Y(),2) );
		double norme_AM = sqrt( pow((A.X()-Point_intersection.X()),2) + pow(A.Y()-Point_intersection.Y(),2) ) ;
		if (norme_AM <= norme_AB){
			math::Vector3d Vect_AB( B.X()-A.X(), B.Y()-A.Y(), 0 );
			math::Vector3d Vect_AM( Point_intersection.X()-A.X(), Point_intersection.Y()-A.Y(), 0 );
			Vect_AB.normalize();
			Vect_AM.normalize();
			if (Vect_AB == Vect_AM){
				M_sur_seg1 = true ;
			}
		}
		// -> Test segment 2 [C,D] si M est bien sur le premier segment
		if (M_sur_seg1){
			double norme_CD = sqrt( pow((C.X()-D.X()),2) + pow(C.Y()-D.Y(),2) );
			double norme_CM = sqrt( pow((C.X()-Point_intersection.X()),2) + pow(C.Y()-Point_intersection.Y(),2) ) ;
			if (norme_CM <= norme_CD){
				math::Vector3d Vect_CD( D.X()-C.X(), D.Y()-C.Y(), 0 );
				math::Vector3d Vect_CM( Point_intersection.X()-C.X(), Point_intersection.Y()-C.Y(), 0 );
				Vect_CD.normalize();
				Vect_CM.normalize();
				if (Vect_CD == Vect_CM){
					M_sur_seg2 = true ;
					intersection = true ;
				}
			}
		}
	}
	return intersection;
}
/*------------------------------------------------------------------------*/