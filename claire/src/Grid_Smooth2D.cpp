//
// Created by claire on 26/11/2021.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/Grid_Smooth2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
Grid_Smooth2D::Grid_Smooth2D(Blocking2D *AMesh,
                   const int ANbIterations) {
	m_mesh = AMesh;
	m_nb_max_iterations = ANbIterations;

}
/*------------------------------------------------------------------------*/
void Grid_Smooth2D::setNbIterations(const int ANbIterations)
{
	m_nb_max_iterations=ANbIterations;

}
/*------------------------------------------------------------------------*/
Grid_Smooth2D::STATUS Grid_Smooth2D::execute()
{

	// Initialization of a pointer to a variable to store the old coordinates of each node
	Variable<math::Point> *old_coords = nullptr;
	old_coords = m_mesh->newVariable<math::Point,GMDS_NODE>("old_coords") ;

	for (int iteration = 1; iteration <= m_nb_max_iterations; iteration++) {

		// Loop on the node ids of the mesh to store the coordinates of each node
		for (auto n_id : m_mesh->nodes()) {
			old_coords->set(n_id, m_mesh->get<Node>(n_id).point());
		}

		// Loop on the blocks of the mesh
		for (auto f_id : m_mesh->faces()) {
			Blocking2D::Block bi = m_mesh->block(f_id) ;
			int Nx = bi.getNbDiscretizationI();
			int Ny = bi.getNbDiscretizationJ();
			std::cout << "Nx = " << Nx << std::endl;
			std::cout << "Ny = " << Ny << std::endl;

			// Loop on the inner nodes of the bi block
			for (int i = 1; i < Nx - 1; i++) {
				for (int j = 1; j < Ny - 1; j++) {

					// Get the stencil points
					math::Point p00 = bi(i - 1, j - 1).point();
					math::Point p01 = bi(i - 1, j).point();
					math::Point p02 = bi(i - 1, j + 1).point();
					math::Point p10 = bi(i, j - 1).point();
					math::Point p11 = bi(i, j).point();
					math::Point p12 = bi(i, j + 1).point();
					math::Point p20 = bi(i + 1, j - 1).point();
					math::Point p21 = bi(i + 1, j).point();
					math::Point p22 = bi(i + 1, j + 1).point();

					// Compute the 6 points for the Yao Smoother
					math::Point V1 = FindMidBranche(p00, p01, p02);
					math::Point V2 = FindMidBranche(p10, p11, p12);
					math::Point V3 = FindMidBranche(p20, p21, p22);

					math::Point H1 = FindMidBranche(p00, p10, p20);
					math::Point H2 = FindMidBranche(p01, p11, p21);
					math::Point H3 = FindMidBranche(p02, p12, p22);

					// Finding the intersection between the 4 segments
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

					if (intersection_trouvee) {
						bi(i, j).setPoint(M);
					}
				}
			}
	}

	}

	return Grid_Smooth2D::SUCCESS;
}

/*------------------------------------------------------------------------*/











/*------------------------------------------------------------------------*/
// Fonction FindMidBranche : Si on considère une branche composée de 3 points, cette fonction retourne le point
//                            positionné au milieu de cette branche
// En entrée : A, B, C -> 3 points
// En sortie : le point milieu
math::Point Grid_Smooth2D::FindMidBranche(const math::Point A, const math::Point B, const math::Point C) {
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

	return Point_Milieu;
}
/*------------------------------------------------------------------------*/
