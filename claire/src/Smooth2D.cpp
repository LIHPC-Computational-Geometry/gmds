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

	// Vérification : est-ce que le maillage est structuré ?
	bool CheckMesh = CheckStructuredMesh();

	if (CheckMesh) {
		// for all free nodes, we build and store stencils locally
		buildStencils();

		for (int iteration = 1; iteration <= m_nb_max_iterations; iteration++) {
			Mesh *new_mesh = m_mesh;
			Node n_new;
			for (auto n_id : m_free_nodes) {
				math::Point H1, H2, H3;
				math::Point V1, V2, V3;
				math::Point A, B, C;
				Node noeud_traite;
				Node noeud_voisin;
				// std::cout << "Noeud :" << n_id << std::endl;
				noeud_traite = m_mesh->get<Node>(n_id);

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

				// Recherche de l'intersection entre les 4 segments
				bool intersection_trouvee(false);
				math::Point M(0, 0, 0);
				math::Segment Seg_Vert_1(V1, V2);
				math::Segment Seg_Vert_2(V2, V3);
				math::Segment Seg_Hori_1(H1, H2);
				math::Segment Seg_Hori_2(H2, H3);

				intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_1, M);
				if (!intersection_trouvee) {
					intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_2, M);
				}
				if (!intersection_trouvee) {
					intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_1, M);
				}
				if (!intersection_trouvee) {
					intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_2, M);
				}

				/*
				// PRINT POUR DEBUG
				if (iteration==m_nb_max_iterations) {
				   std::cout << "--------------------------------" << std::endl;
				   std::cout << "Iterations :" << iteration << std::endl;
				   std::cout << "Noeud :" << n_id << std::endl;
				   std::cout << "Intersection trouvee ?" << intersection_trouvee << std::endl;
				   std::cout << "Pos actuelle :" << noeud_traite.point() << std::endl;
				   std::cout << "Nouvelle position :" << M << std::endl;
				}
				// FIN DES PRINT POUR DEBUG
				*/

				n_new = new_mesh->get<Node>(n_id);
				if (intersection_trouvee) {
					n_new.setPoint(M);
				}
			}

			m_mesh = new_mesh;
		}
	}

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
// ------------ PERTURBATION ALEATOIRE DU MAILLAGE ------------
void Smooth2D::PerturbationMaillage(const Variable<int>* var_bnd, const double dx, const double dy) {
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
	math::Vector3d Vec_AB = B-A ;
	math::Vector3d Vec_BC = C-B ;
	double norme_1 = Vec_AB.norm() ;
	double norme_2 = Vec_BC.norm() ;
	double norme_branche = norme_1 + norme_2 ;
	double norme_milieu = norme_branche / 2.0 ;

	if (norme_milieu <= norme_1){
		Vec_AB.normalize();
		Point_Milieu = A + norme_milieu*Vec_AB ;
	}
	else if (norme_milieu > norme_1){
		math::Vector3d Vec_CB = - Vec_BC ;
		Vec_CB.normalize();
		Point_Milieu = C + norme_milieu*Vec_CB ;
	}

	/*
	std::cout << "Point milieu :" << Point_Milieu << std::endl ;
	std::cout << "--------------------------------------" << std::endl;
	 */
	return Point_Milieu;
}
/*------------------------------------------------------------------------*/










/*------------------------------------------------------------------------*/
// Fonction CheckStructuredMesh : Prend un maillage en entrée et vérifie si il est structuré ou non
// En entrée :
// En sortie :
bool Smooth2D::CheckStructuredMesh() {
	bool checkmesh(true);
	// Vérifie que chaque noeud intérieur a bien 4 faces adjacentes
	for(auto n_id:m_free_nodes) {
		if (checkmesh) {
			Node current_node = m_mesh->get<Node>(n_id);
			std::vector<Face> current_faces = current_node.get<Face>();
			std::cout << "Current faces :" << current_faces[0] << std::endl;
			int Nbr_Faces = current_faces.size();
			if (Nbr_Faces != 4) {
				checkmesh = false;
				std::cout << "Numéro du noeud : " << n_id << std::endl;
				std::cout << "Nbr de faces : " << Nbr_Faces << std::endl;
				std::cout << "-----------" << std::endl;
			}
		}
	}

	// Vérifie que chaque face est un quad
	for(auto n_id:m_mesh->faces()){
		if (checkmesh) {
			std::cout << "Numéro de la face : " << n_id << std::endl;
			Face current_face = m_mesh->get<Face>(n_id);
			std::vector<Node> current_nodes = current_face.get<Node>();
			std::cout << "Current nodes :" << current_nodes[0] << std::endl;
			int Nbr_Noeuds = current_nodes.size();
			if (Nbr_Noeuds != 4) {
				checkmesh = false;
				std::cout << "Nbr de noeuds : " << Nbr_Noeuds << std::endl;
				std::cout << "-----------" << std::endl;
			}
		}
	}

	std::cout << "Maillage structuré ? " << checkmesh << std::endl ;
	return checkmesh;
}
/*------------------------------------------------------------------------*/
