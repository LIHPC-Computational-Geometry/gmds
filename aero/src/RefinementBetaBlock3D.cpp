//
// Created by rochec on 29/05/24.
//

/*------------------------------------------------------------------------*/
#include <gmds/aero/RefinementBetaBlock3D.h>
#include <gmds/aero/RefinementBeta.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

RefinementBetaBlock3D::RefinementBetaBlock3D(Blocking3D::Block* ABlock,
                                             Blocking3D::BlockFace* Abf,
                                             double Asize_first_edge)
{
	m_Block = ABlock;
	m_bf = Abf;
	m_size_first_edge = Asize_first_edge;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
RefinementBetaBlock3D::STATUS
RefinementBetaBlock3D::execute()
{
	if (m_Block->isFaceI0(m_bf->id()))
	{
		for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
		{
			for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
			{
				std::vector<math::Point> pos;
				for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
				{
					pos.push_back((*m_Block)(i,j,k).point());
				}
				RefinementBeta algo_ref = RefinementBeta(pos, m_size_first_edge);
				algo_ref.execute();
				std::vector<math::Point> new_pos = algo_ref.GetNewPositions();
				for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
				{
					(*m_Block)(i,j,k).setPoint(new_pos[i]);
				}
			}
		}
	}
	else if (m_Block->isFaceImax(m_bf->id()))
	{
		for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
		{
			for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
			{
				std::vector<math::Point> pos;
				for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
				{
					pos.push_back((*m_Block)(m_Block->getNbDiscretizationI()-(i+1),j,k).point());
				}
				RefinementBeta algo_ref = RefinementBeta(pos, m_size_first_edge);
				algo_ref.execute();
				std::vector<math::Point> new_pos = algo_ref.GetNewPositions();
				for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
				{
					(*m_Block)(m_Block->getNbDiscretizationI()-(i+1),j,k).setPoint(new_pos[i]);
				}
			}
		}
	}
	else if (m_Block->isFaceJ0(m_bf->id()))
	{
		for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
		{
			for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
			{
				std::vector<math::Point> pos;
				for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
				{
					pos.push_back((*m_Block)(i,j,k).point());
				}
				RefinementBeta algo_ref = RefinementBeta(pos, m_size_first_edge);
				algo_ref.execute();
				std::vector<math::Point> new_pos = algo_ref.GetNewPositions();
				for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
				{
					(*m_Block)(i,j,k).setPoint(new_pos[j]);
				}
			}
		}
	}
	else if (m_Block->isFaceJmax(m_bf->id()))
	{
		for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
		{
			for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
			{
				std::vector<math::Point> pos;
				for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
				{
					pos.push_back((*m_Block)(i,m_Block->getNbDiscretizationJ()-(j+1),k).point());
				}
				RefinementBeta algo_ref = RefinementBeta(pos, m_size_first_edge);
				algo_ref.execute();
				std::vector<math::Point> new_pos = algo_ref.GetNewPositions();
				for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
				{
					(*m_Block)(i,m_Block->getNbDiscretizationJ()-(j+1),k).setPoint(new_pos[j]);
				}
			}
		}
	}
	else if (m_Block->isFaceK0(m_bf->id()))
	{
		for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
		{
			for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
			{
				std::vector<math::Point> pos;
				for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
				{
					pos.push_back((*m_Block)(i,j,k).point());
				}
				RefinementBeta algo_ref = RefinementBeta(pos, m_size_first_edge);
				algo_ref.execute();
				std::vector<math::Point> new_pos = algo_ref.GetNewPositions();
				for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
				{
					(*m_Block)(i,j,k).setPoint(new_pos[k]);
				}
			}
		}
	}
	else if (m_Block->isFaceKmax(m_bf->id()))
	{
		for (auto i=0;i<m_Block->getNbDiscretizationI();i++)
		{
			for (auto j=0;j<m_Block->getNbDiscretizationJ();j++)
			{
				std::vector<math::Point> pos;
				for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
				{
					pos.push_back((*m_Block)(i,j,m_Block->getNbDiscretizationK()-(k+1)).point());
				}
				RefinementBeta algo_ref = RefinementBeta(pos, m_size_first_edge);
				algo_ref.execute();
				std::vector<math::Point> new_pos = algo_ref.GetNewPositions();
				for (auto k=0;k<m_Block->getNbDiscretizationK();k++)
				{
					(*m_Block)(i,j,m_Block->getNbDiscretizationK()-(k+1)).setPoint(new_pos[k]);
				}
			}
		}
	}

	return RefinementBetaBlock3D::SUCCESS;
}
/*------------------------------------------------------------------------*/