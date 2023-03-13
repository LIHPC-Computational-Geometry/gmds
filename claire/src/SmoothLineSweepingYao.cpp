//
// Created by rochec on 21/07/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/SmoothLineSweepingYao.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

SmoothLineSweepingYao::SmoothLineSweepingYao(Blocking2D::Block *AB, int Anb_max_it, double Atheta) :
  AbstractSmoothLineSweeping_2D(AB, Anb_max_it, Atheta)
{

}
/*------------------------------------------------------------------------*/
math::Point SmoothLineSweepingYao::ComputeNewPosition(int i, int j)
{
	// Compute the 6 points for the Yao Smoother
	math::Point V1 = math::Utils::WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i,j-1).point(), (*m_B)(i+1,j-1).point(), 0.5);
	math::Point V2 = math::Utils::WeightedPointOnBranch((*m_B)(i-1,j).point(), (*m_B)(i,j).point(), (*m_B)(i+1,j).point(), 0.5);
	math::Point V3 = math::Utils::WeightedPointOnBranch((*m_B)(i-1,j+1).point(), (*m_B)(i,j+1).point(), (*m_B)(i+1,j+1).point(), 0.5);

	math::Point H1 = math::Utils::WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i-1,j).point(), (*m_B)(i-1,j+1).point(), 0.5);
	math::Point H2 = math::Utils::WeightedPointOnBranch((*m_B)(i,j-1).point(), (*m_B)(i,j).point(), (*m_B)(i,j+1).point(), 0.5);
	math::Point H3 = math::Utils::WeightedPointOnBranch((*m_B)(i+1,j-1).point(), (*m_B)(i+1,j).point(), (*m_B)(i+1,j+1).point(), 0.5);

	// Finding the intersection between the 4 segments
	math::Point M;
	math::Segment Seg_Vert_1(V1, V2);
	math::Segment Seg_Vert_2(V2, V3);
	math::Segment Seg_Hori_1(H1, H2);
	math::Segment Seg_Hori_2(H2, H3);

	bool intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_1, M);
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_2, M);
	}
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_1, M);
	}
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_2, M);
	}

	math::Point P_ideal = (*m_B)(i,j).point() ;
	if (intersection_trouvee) {
		P_ideal = M;
	}

	return P_ideal;
}
/*------------------------------------------------------------------------*/