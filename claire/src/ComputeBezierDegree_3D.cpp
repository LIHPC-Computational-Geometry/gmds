//
// Created by rochec on 17/01/24.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/ComputeBezierDegree_3D.h>
#include <gmds/math/BezierSurface.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/claire/Utils.h>
#include <gmds/ig/Mesh.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

ComputeBezierDegree_3D::ComputeBezierDegree_3D(Mesh *AMeshHex,
                                               cad::FACManager* AManager,
                                               cad::GeomMeshLinker* ALinker_HG,
                                               std::vector<TCellID>* Afaces_to_fit,
                                               double Amax_error,
                                               int Amax_degree) :
  m_meshHex(AMeshHex),
  m_manager(AManager),
  m_linker_HG(ALinker_HG),
  m_faces_to_fit(Afaces_to_fit),
  m_max_error(Amax_error),
  m_max_degree(Amax_degree),
  m_final_degree(1)
{

}


/*------------------------------------------------------------------------*/
ComputeBezierDegree_3D::STATUS
ComputeBezierDegree_3D::execute()
{
	int final_degree(1);
	int iter(0);
	//std::cout << "Nbr faces : " << m_faces_to_fit->size() << std::endl;
	while( (iter < m_faces_to_fit->size()) && (final_degree < m_max_degree) )
	{
		//std::cout << "Face " << iter << std::endl;
		TCellID f_id = (*m_faces_to_fit)[iter];
		Face f = m_meshHex->get<Face>(f_id);
		cad::GeomSurface* surf = m_manager->getSurface( m_linker_HG->getGeomId(f) ) ;
		int loc_degree(final_degree);

		Array2D<math::Point> init_ctrl_pts = computeBezierCtrlPtstoInterpolateSurface(f, loc_degree, surf);
		math::BezierSurface init_bs(init_ctrl_pts);
		double loc_error = computeErrorBtwBezierSurfaceandGeomSurface(init_bs, surf);

		while ( (loc_degree<m_max_degree) && (loc_error>m_max_error) )
		{
			Array2D<math::Point> ctrl_pts = computeBezierCtrlPtstoInterpolateSurface(f, loc_degree, surf);
			math::BezierSurface bs(ctrl_pts);
			loc_error = computeErrorBtwBezierSurfaceandGeomSurface(bs, surf);
			//std::cout << "Error " << loc_error << std::endl;
			if (loc_error>m_max_error)
			{
				loc_degree++;
			}
		}
		final_degree = loc_degree;

		iter++;

	}

	m_final_degree = final_degree;

	return ComputeBezierDegree_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
int
ComputeBezierDegree_3D::getDegree()
{
	return m_final_degree;
}
/*------------------------------------------------------------------------*/
Array2D<math::Point>
ComputeBezierDegree_3D::computeBezierCtrlPtstoInterpolateSurface(Face f, int degree, cad::GeomSurface* surf)
{
	Array2D<math::Point> ctrl_pts(degree+1,degree+1);
	Array2D<math::Point> interpolate_pts(degree+1,degree+1);

	std::vector<Node> f_nodes = f.get<Node>();
	if (f_nodes.size() != 4)
	{
		std::cout << "WARNING in ComputeBezierDegree_3D: a face is not a quad." << std::endl;
		std::cout << "Face " << f.id() << std::endl;
	}

	// Init the array of points to interpolate using the nodes of the face and transfinite interpolation
	interpolate_pts(0,0) = f_nodes[0].point() ;
	interpolate_pts(degree,0) = f_nodes[1].point() ;
	interpolate_pts(degree,degree) = f_nodes[2].point() ;
	interpolate_pts(0,degree) = f_nodes[3].point() ;

	for (int i=1;i<degree;i++)
	{
		interpolate_pts(i,0) = interpolate_pts(0,0) + i*(interpolate_pts(degree,0)- interpolate_pts(0,0))/(float(degree)) ;
		interpolate_pts(i,degree) = interpolate_pts(0,degree) + i*(interpolate_pts(degree,degree)- interpolate_pts(0,degree))/(float(degree)) ;
	}
	for (int j=1;j<degree;j++)
	{
		interpolate_pts(0,j) = interpolate_pts(0,0) + j*(interpolate_pts(0,degree)- interpolate_pts(0,0))/(float(degree)) ;
		interpolate_pts(degree,j) = interpolate_pts(degree,0) + j*(interpolate_pts(degree,degree)- interpolate_pts(degree,0))/(float(degree)) ;
	}

	math::TransfiniteInterpolation::computeQuad(interpolate_pts);

	// Project the nodes onto the surface to interpolate
	for (int i=0;i<=degree;i++)
	{
		for (int j=0;j<=degree;j++)
		{
			math::Point p = interpolate_pts(i,j);
			surf->project(p);
			interpolate_pts(i,j) = p;
		}
	}

	// Now we compute the control points to interpolate the array interpolate_pts
	Eigen::MatrixXd mat_B((degree+1)*(degree+1), (degree+1)*(degree+1));
	Eigen::VectorXd ctrl_points_x((degree+1)*(degree+1));
	Eigen::VectorXd ctrl_points_y((degree+1)*(degree+1));
	Eigen::VectorXd ctrl_points_z((degree+1)*(degree+1));
	Eigen::VectorXd interp_points_x((degree+1)*(degree+1));
	Eigen::VectorXd interp_points_y((degree+1)*(degree+1));
	Eigen::VectorXd interp_points_z((degree+1)*(degree+1));
	for (int i=0;i<=degree;i++)
	{
		for (int j=0;j<=degree;j++)
		{
			interp_points_x[i+j*(degree+1)] = interpolate_pts(i,j).X();
			interp_points_y[i+j*(degree+1)] = interpolate_pts(i,j).Y();
			interp_points_z[i+j*(degree+1)] = interpolate_pts(i,j).Z();
		}
	}

	// Matrix Assembly
	for (int i=0;i<=degree;i++)
	{
		for (int j=0;j<=degree;j++)
		{
			double u = 1.0*i/degree ;
			double v = 1.0*j/degree ;
			for (int i2=0;i2<=degree;i2++)
			{
				for (int j2=0;j2<=degree;j2++)
				{
					double Bi2n = math::Utils::BernsteinPolynomial(degree, i2, u) ;
					double Bj2m = math::Utils::BernsteinPolynomial(degree, j2, v) ;
					mat_B(i+j*(degree+1), i2+j2*(degree+1)) = Bi2n*Bj2m ;
				}
			}
		}
	}

	Eigen::MatrixXd mat_B_inv = mat_B.inverse();
	ctrl_points_x = mat_B_inv*interp_points_x;
	ctrl_points_y = mat_B_inv*interp_points_y;
	ctrl_points_z = mat_B_inv*interp_points_z;

	for (int i=0;i<=degree;i++)
	{
		for (int j=0;j<=degree;j++)
		{
			math::Point p;
			p.setX( ctrl_points_x(i+j*(degree+1)) );
			p.setY( ctrl_points_y(i+j*(degree+1)) );
			p.setZ( ctrl_points_z(i+j*(degree+1)) );
			ctrl_pts(i,j) = p;
		}
	}

	return ctrl_pts;

}
/*------------------------------------------------------------------------*/
double
ComputeBezierDegree_3D::computeErrorBtwBezierSurfaceandGeomSurface(math::BezierSurface Abs, cad::GeomSurface* surf)
{
	int sample(20);
	double err(0.0);

	for (int i=1;i<sample;i++)
	{
		for (int j=1;j<sample;j++)
		{
			double u(float(i)/float(sample));
			double v(float(j)/float(sample));
			math::Point p = Abs(u,v);
			surf->project(p);
			double loc_err = (Abs(u,v)-p).norm() ;
			if (loc_err > err)
			{
				err = loc_err;
				//std::cout << "u,v " << u << ", " << v << ": " << loc_err << std::endl;
			}
		}
	}

	// Check ctrl pts on edges
	//std::cout << "check on edges" << std::endl;
	for (int i=1;i<sample;i++)
	{
		double u(float(i)/float(sample));
		math::Point p = Abs(u,0);
		surf->project(p);
		double loc_err = (Abs(u,0)-p).norm() ;
		if (loc_err > err)
		{
			err = loc_err;
			//std::cout << "u,v " << u << ", " << 0 << ": " << loc_err << std::endl;
		}
		p = Abs(u,1);
		surf->project(p);
		loc_err = (Abs(u,1)-p).norm() ;
		if (loc_err > err)
		{
			err = loc_err;
			//std::cout << "u,v " << u << ", " << 1 << ": " << loc_err << std::endl;
		}
	}

	for (int j=1;j<sample;j++)
	{
		double v(float(j)/float(sample));
		math::Point p = Abs(0,v);
		surf->project(p);
		double loc_err = (Abs(0,v)-p).norm() ;
		if (loc_err > err)
		{
			err = loc_err;
			//std::cout << "u,v " << 0 << ", " << v << ": " << loc_err << std::endl;
		}
		p = Abs(1,v);
		surf->project(p);
		loc_err = (Abs(1,v)-p).norm() ;
		if (loc_err > err)
		{
			err = loc_err;
			//std::cout << "u,v " << 1 << ", " << v << ": " << loc_err << std::endl;
		}
	}

	//std::cout << "Max error on face : " << err << std::endl;
	return err;

}
/*------------------------------------------------------------------------*/