//
// Created by rochec on 28/11/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/ControlPointsSmoothing_2D.h>
#include <gmds/claire/Utils.h>
#include <gmds/math/BezierCurve.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
ControlPointsSmoothing_2D::ControlPointsSmoothing_2D(Blocking2D* ACtrlPoints) :
m_CtrlPts(ACtrlPoints)
{

}
/*------------------------------------------------------------------------*/
ControlPointsSmoothing_2D::STATUS ControlPointsSmoothing_2D::execute()
{
	for (auto bloc:m_CtrlPts->allBlocks())
	{
		std::vector<math::Point> Pts;

		//First edge: j=0
		for (int i=0;i<bloc.getNbDiscretizationI();i++)
		{
			Pts.push_back(bloc(i,0).point());
		}

		for (int i = 2; i < bloc.getNbDiscretizationI()-2; i++)
		{
			math::Vector3d tangent = math::Utils::DerivativeBezierCurve(Pts, 1.0*i/(bloc.getNbDiscretizationI()-1.0)) ;
			math::Vector3d normal ;
			normal.setX(-tangent.Y());
			normal.setY(tangent.X());
			normal.normalize();

			math::Vector3d v = bloc(i,1).point()-bloc(i,0).point() ;

			if (v.dot(normal) >= 0)
			{
				math::Point p = bloc(i,0).point() + v.norm()*normal ;
				bloc(i,1).setPoint(p);
			}
			else
			{
				math::Point p = bloc(i,0).point() + (-v.norm())*normal ;
				bloc(i,1).setPoint(p);
			}
		}

		//Second edge: j=bloc.getNbDiscretizationJ()-1
		Pts.clear();
		for (int i=0;i<bloc.getNbDiscretizationI();i++)
		{
			Pts.push_back(bloc(i,bloc.getNbDiscretizationJ()-1).point());
		}

		for (int i = 2; i < bloc.getNbDiscretizationI()-2; i++)
		{
			math::Vector3d tangent = math::Utils::DerivativeBezierCurve(Pts, 1.0*i/(bloc.getNbDiscretizationI()-1.0)) ;
			math::Vector3d normal ;
			normal.setX(-tangent.Y());
			normal.setY(tangent.X());
			normal.normalize();

			math::Vector3d v = bloc(i,bloc.getNbDiscretizationJ()-2).point()-bloc(i,bloc.getNbDiscretizationJ()-1).point() ;

			if (v.dot(normal) >= 0)
			{
				math::Point p = bloc(i,bloc.getNbDiscretizationJ()-1).point() + v.norm()*normal ;
				bloc(i,bloc.getNbDiscretizationJ()-2).setPoint(p);
			}
			else
			{
				math::Point p = bloc(i,bloc.getNbDiscretizationJ()-1).point() + (-v.norm())*normal ;
				bloc(i,bloc.getNbDiscretizationJ()-2).setPoint(p);
			}
		}


		//Third edge: i=0
		Pts.clear();
		for (int j=0;j<bloc.getNbDiscretizationJ();j++)
		{
			Pts.push_back(bloc(0,j).point());
		}

		for (int j = 2; j < bloc.getNbDiscretizationJ()-2; j++)
		{
			math::Vector3d tangent = math::Utils::DerivativeBezierCurve(Pts, 1.0*j/(bloc.getNbDiscretizationJ()-1.0)) ;
			math::Vector3d normal ;
			normal.setX(-tangent.Y());
			normal.setY(tangent.X());
			normal.normalize();

			math::Vector3d v = bloc(1,j).point()-bloc(0,j).point() ;

			if (v.dot(normal) >= 0)
			{
				math::Point p = bloc(0,j).point() + v.norm()*normal ;
				bloc(1,j).setPoint(p);
			}
			else
			{
				math::Point p = bloc(0,j).point() + (-v.norm())*normal ;
				bloc(1,j).setPoint(p);
			}
		}


		//Last edge: i=bloc.getNbDiscretizationI()-1
		Pts.clear();
		for (int j=0;j<bloc.getNbDiscretizationJ();j++)
		{
			Pts.push_back(bloc(bloc.getNbDiscretizationI()-1,j).point());
		}

		for (int j = 2; j < bloc.getNbDiscretizationJ()-2; j++)
		{
			math::Vector3d tangent = math::Utils::DerivativeBezierCurve(Pts, 1.0*j/(bloc.getNbDiscretizationJ()-1.0)) ;
			math::Vector3d normal ;
			normal.setX(-tangent.Y());
			normal.setY(tangent.X());
			normal.normalize();

			math::Vector3d v = bloc(bloc.getNbDiscretizationI()-2,j).point()-bloc(bloc.getNbDiscretizationI()-1,j).point() ;

			if (v.dot(normal) >= 0)
			{
				math::Point p = bloc(bloc.getNbDiscretizationI()-1,j).point() + v.norm()*normal ;
				bloc(bloc.getNbDiscretizationI()-2,j).setPoint(p);
			}
			else
			{
				math::Point p = bloc(bloc.getNbDiscretizationI()-1,j).point() + (-v.norm())*normal ;
				bloc(bloc.getNbDiscretizationI()-2,j).setPoint(p);
			}
		}


	}

	return ControlPointsSmoothing_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/