/*----------------------------------------------------------------------------*/
#include "gmds/smoothy/EllipticSmoothing.h"
#include <gmds/smoothy/HLBFGSWrapper.h>
/*----------------------------------------------------------------------------*/
#include <cmath>
#include <iostream>
#include<limits>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
using namespace gmds;
/*----------------------------------------------------------------------------*/
double chi(double eps, double det) {
	const double eps2 = eps * eps;
    if (det > 0)
        return (det + std::sqrt(eps2 + det * det)) * .5;
    return .5 * eps2 / (std::sqrt(eps2 + det * det) - det);
}
/*----------------------------------------------------------------------------*/
double chi_deriv(double eps, double det) {
    return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}
/*----------------------------------------------------------------------------*/
void display_options(const EllipticSmoothingOptions & opt) {
    std::cerr << "Options: " << std::endl;
    std::cerr << "-theta = " << opt.theta << std::endl;
    std::cerr << "-maxiter = " << opt.maxiter << std::endl;
    std::cerr << "-bfgs_threshold = " << opt.bfgs_threshold << std::endl;
    std::cerr << "-bfgs_maxiter = " << opt.bfgs_maxiter << std::endl;
    std::cerr << "-debug = " << opt.debug << std::endl;
    std::cerr << "-eps_from_theorem = " << opt.eps_from_theorem << std::endl;
    std::cerr << "-stopping_when_static = " << opt.stopping_when_static << std::endl;
    std::cerr << "-static_threshold = " << opt.static_threshold << std::endl;
}
/*----------------------------------------------------------------------------*/
math::Matrix33 dual_basis(const math::Matrix33& J) {
	return
	   {
	      J[1].Y() * J[2].Z() - J[1].Z() * J[2].Y(),
	      J[1].Z() * J[2].X() - J[1].X() * J[2].Z(),
	      J[1].X() * J[2].Y() - J[1].Y() * J[2].X(),
	      J[0].Z() * J[2].Y() - J[0].Y() * J[2].Z(),
	      J[0].X() * J[2].Z() - J[0].Z() * J[2].X(),
	      J[0].Y() * J[2].X() - J[0].X() * J[2].Y(),
	      J[0].Y() * J[1].Z() - J[0].Z() * J[1].Y(),
	      J[0].Z() * J[1].X() - J[0].X() * J[1].Z(),
	      J[0].X() * J[1].Y() - J[0].Y() * J[1].X()
	   };
}
/*----------------------------------------------------------------------------*/
void EllipticSmoothing2D::set_triangle_ref(const int t, const std::array<double, 6> &ref) {
    assert(t < m_nb_tri);
    mat22 M = { ref[2] - ref[0], ref[4] - ref[0], ref[3] - ref[1], ref[5] - ref[1] };
    double detM = (M[0] * M[3] - M[1] * M[2]);
	 m_area[t] = 0.5 * detM;
    double invdetM = 1. / detM;
    mat22 invM = { invdetM * M[3], -invdetM * M[1], -invdetM * M[2], invdetM * M[0] };
    refs_grad_[t] = { -invM[0] - invM[2], -invM[1] - invM[3], invM[0], invM[1], invM[2], invM[3] };
}
/*----------------------------------------------------------------------------*/
void EllipticSmoothing2D::set_triangle_ref(const std::vector<std::array<double, 6>> &refs) {
    assert(refs.size() == m_nb_tri);
    for(auto i_tri=0;i_tri<m_nb_tri;i_tri++)
		 set_triangle_ref(i_tri, refs[i_tri]);
}
/*----------------------------------------------------------------------------*/
void EllipticSmoothing2D::set_area(const int t, const double area) {
    assert(t < m_nb_tri);
	 m_area[t] = area;
}
/*----------------------------------------------------------------------------*/
void EllipticSmoothing2D::init() {
    refs_grad_.resize(m_nb_tri);
	 m_area.resize(m_nb_tri);
	 for(auto i_tri=0;i_tri<m_nb_tri;i_tri++)  {
        refs_grad_[i_tri] = { { 0,-1, std::sqrt(3.) / 2.,.5, -std::sqrt(3.) / 2.,.5 } };
		  m_area[i_tri] = std::sqrt(3) / 4;
    }
	 m_J.resize(m_nb_tri);
	 m_K.resize(m_nb_tri);
	 m_det.resize(m_nb_tri);
}
/*----------------------------------------------------------------------------*/
void EllipticSmoothing2D::evaluate_jacobian() {
    if (m_options.debug > 3) std::cerr << "evaluate the jacobian...";
	 m_detmin = std::numeric_limits<double>::max();
	 m_nb_inverted = 0;

	 for(auto i_t=0;i_t<m_nb_tri;i_t++){
		 m_J[i_t] = { 0,0,0,0 };
		 for(auto i=0;i<3;i++) {
			 for(auto d=0;d<2;d++){
				 const double v = var_.value(6 * i_t + 2 * i + d);
				 for(auto dj=0;dj<2;dj++)
					 m_J[i_t][2 * d + dj] += refs_grad_[i_t][2 * i + dj] * v;
			 }
		 }

		m_det[i_t] = m_J[i_t][0] * m_J[i_t][3] - m_J[i_t][1] * m_J[i_t][2];
		if(m_det[i_t] < m_detmin)
			m_detmin = m_det[i_t];
		if(m_det[i_t] <= 0.) ++m_nb_inverted;

		m_K[i_t] = { {m_J[i_t][3], -m_J[i_t][2], -m_J[i_t][1], m_J[i_t][0] } };// dual basis
    }
    if (m_options.debug > 3) std::cerr << "ok" << std::endl;
}
/*----------------------------------------------------------------------------*/
double EllipticSmoothing2D::evaluate_energy() {
    evaluate_jacobian();
    double E = 0;
    for(auto t=0;t<m_nb_tri;t++) {
        double chi_ = chi(m_eps, m_det[t]);
        double f = 0;
        for(auto i=0;i<4;i++)
			  f += m_J[t][i] * m_J[t][i];
        double g = 1 + m_det[t] * m_det[t];
        E += ((1. - m_options.theta) * f + m_options.theta * g) * m_area[t] / chi_;
    }
    return E;
}
/*----------------------------------------------------------------------------*/
void EllipticSmoothing2D::standard_elliptic_energy(const std::vector<double>& X, double& F, std::vector<double>& G) {
    evaluate_jacobian();

	  for(auto t=0;t<m_nb_tri;t++) {
        double c1 = chi(m_eps, m_det[t]);
        double c2 = chi_deriv(m_eps, m_det[t]);

        double f = 0;
		  for(auto i=0;i<4;i++)
			  f += m_J[t][i] * m_J[t][i];
        f /= c1;
        double g = (1 + m_det[t] * m_det[t]) / c1;
        F += ((1. - m_options.theta) * f + m_options.theta * g) * m_area[t];

        const double df_mul_a = 2. / c1;
        const double df_mul_b = f * c2 / c1;
        const double dg_mul = (2 * m_det[t] - g * c2) / c1;

        for(auto d=0;d<2;d++) {
            double a1 = m_J[t][2 * d], a2 = m_J[t][2 * d + 1]; // tangent basis
            double b1 = m_K[t][2 * d], b2 = m_K[t][2 * d + 1]; // dual basis
            double dfda1 = a1 * df_mul_a - b1 * df_mul_b, dfda2 = a2 * df_mul_a - b2 * df_mul_b;
            double dgda1 = b1 * dg_mul, dgda2 = b2 * dg_mul;
            for(auto i =0;i<3;i++) {
                double gradi = (dfda1 * (1. - m_options.theta) + dgda1 * m_options.theta) * refs_grad_[t][2 * i]
                             + (dfda2 * (1. - m_options.theta) + dgda2 * m_options.theta) * refs_grad_[t][2 * i + 1];
                gradi *= m_area[t];
                var_.add2grad(6 * t + 2 * i + d, gradi, G);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void EllipticSmoothing2D::standard_elliptic_energy_w_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<EllipticSmoothingMeshGlue::sparse_term>& H) {
	standard_elliptic_energy(X, F, G);
}
/*----------------------------------------------------------------------------*/
bool EllipticSmoothing2D::run_lbfgs(std::vector<double>& X) {
	bool first = true;
	double F0, F1;
		const HLBFGSWrapper::gradient_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
        std::fill(G.begin(), G.end(), 0);
        F = 0;
		  var_.set_reduced_values(X);
        energy(X, F, G);
		if(first) {
			F0 = F1 = F;
			first = false;
		} else F1 = std::min(F1, F);
    };
	   HLBFGSWrapper opt = { func };
		opt.gtol = m_options.bfgs_threshold;
    	//opt.maxiter = options_.bfgs_maxiter;
    	//opt.verbose = options_.debug;;
  //  opt.invH.history_depth = 20;

    //if(!options_.eps_from_theorem && detmin_ > 0.) {
    //    opt.gtol = std::min(options_.bfgs_threshold, .1 * options_.static_threshold);
    //    opt.maxiter = options_.bfgs_maxiter * std::max(1, options_.maxiter / 2);
    //} else {
        opt.gtol = m_options.bfgs_threshold;
        opt.maxiter = m_options.bfgs_maxiter;
    //}
    opt.verbose = m_options.debug > 0;


    opt.run(X);
    if (m_options.debug > 0) std::cerr << " F : " << F0 << "  --->  " << F1 << "    ||   ";
    return true;// (F0 - F1) / F1 > 2.e-7;
}
/*----------------------------------------------------------------------------*/
bool EllipticSmoothing2D::run_newton(std::vector<double>& X) {
    return true;
}
/*----------------------------------------------------------------------------*/
bool EllipticSmoothing2D::execute() {
    if (m_options.debug > 0) {
        std::cerr << "==== Running Elliptic smoother 2d. ====" << std::endl;
        display_options(m_options);
    }
    if (var_.get_nb_reduced_variables() == 0) {
        std::cerr << "No variables to optimize" << std::endl;
        return true;
    }
    double total_area_ = 0;
    for(auto i=0;i<m_nb_tri;i++)
		 total_area_ += m_area[i];
	 for(auto i=0;i<m_nb_tri;i++)
		 m_area[i] /= total_area_;

    evaluate_jacobian();
    const double e0 = 1e-3;
    if (m_options.eps_from_theorem)
		 m_eps = start_eps;
    else if (m_options.barrier)
		 m_eps = 1e-7;
    else
		 m_eps = m_detmin > 0 ? .5*e0 : std::sqrt(e0*e0 + 0.04 * m_detmin * m_detmin);
    std::vector<double> X(var_.get_nb_reduced_variables());
	 var_.get_reduced_values(X);
    bool nullStep = false;

    for(auto iter=0; iter<m_options.maxiter;iter++) {
        const double E_prev = evaluate_energy();
        const double detmin_prev = m_detmin;
        if (m_options.debug > 0) std::cerr << "iteration #" << iter << ":    eps: " << m_eps << " detmin: " << m_detmin << " ninv: " << m_nb_inverted << std::endl;

		bool better;
        if (m_options.use_newton) better = run_newton(X);
        else better = run_lbfgs(X);
		  var_.set_reduced_values(X);

        const double E = evaluate_energy();
        if (m_options.debug > 0) std::cerr << " E : " << E_prev << "  --->  " << E << std::endl;
        if (m_options.eps_from_theorem) {
            const double sigma = std::max(1. - E / E_prev, 1e-1);
            if (m_detmin >= 0)
				   m_eps *= (1 - sigma);
            else {
				const double det_eps_norm = std::sqrt(m_detmin * m_detmin + m_eps * m_eps);
				m_eps *= 1 - (sigma * det_eps_norm) / (std::abs(m_detmin) + det_eps_norm);
			}
        } else if (detmin_prev > 0. && m_detmin > 0.) {
            std::cerr << "Stopping as detmin > 0 while not using the eps from theorem" << std::endl;
            break;
        } else if(!m_options.barrier) {
			  m_eps = std::min(.995* m_eps, m_detmin > 0 ? .5*e0 : std::sqrt(e0*e0 + 0.04 * m_detmin * m_detmin));
        }
        if (m_detmin > 0 || m_options.stopping_when_static) {
            if ((E_prev - E) / E < m_options.static_threshold) break;
        } else if (!better) {
            if(nullStep) break;
            nullStep = true;
        } else nullStep = false;
    }
    if (m_options.debug > 0) std::cerr << "E: " << evaluate_energy() << " detmin: " << m_detmin << " ninv: " << m_nb_inverted << std::endl;
    return !m_nb_inverted;
}
/*----------------------------------------------------------------------------*/
void Elliptic_smoother_3D::set_tet_ref(const int t, const std::array < math::Vector3d , 4> tet_ref) {
	assert(t < N_tets_);
	math::Matrix33   M = { tet_ref[1] - tet_ref[0], tet_ref[2] - tet_ref[0], tet_ref[3] - tet_ref[0] };
	math::Matrix33 invM = M.inverse();
	double detM = M.det();
	vol_[t] = 1./6. * detM;
	invM = invM.transpose();
	refs_grad_[t] = { -invM[0] - invM[1] - invM[2], invM[0], invM[1], invM[2] };
}
/*----------------------------------------------------------------------------*/
void Elliptic_smoother_3D::set_tet_ref(const std::vector<std::array < math::Vector3d , 4> > tets_ref) {
	assert(tets_ref.size() == N_tets_);
	for (int t = 0; t < N_tets_; t++) set_tet_ref(t, tets_ref[t]);
}
/*----------------------------------------------------------------------------*/
void Elliptic_smoother_3D::init() {
    refs_grad_.resize(N_tets_);
    vol_.resize(N_tets_);
    for(auto t=0;t<N_tets_;t++) {
        constexpr double a = 0.70710678118; // invsqrt of 2, not possible in constexpr...
		  refs_grad_[t] = { math::Vector3d({-a,-a,-a}), math::Vector3d({a,a,-a}), math::Vector3d({-a,a,a}), math::Vector3d({a,-a,a}) };
        vol_[t] = 1./6.*a;
    }
    J_.resize(N_tets_);
    K_.resize(N_tets_);
    det_.resize(N_tets_);
}
/*----------------------------------------------------------------------------*/
void Elliptic_smoother_3D::evaluate_jacobian() {
    if (options_.debug > 3) std::cerr << "evaluate the jacobian...";
    detmin_ = std::numeric_limits<double>::max();
    ninverted_ = 0;
    for(auto t=0;t<N_tets_;t++) {
		 J_[t] = { math::Vector3d({0,0,0}), math::Vector3d({0,0,0}),math::Vector3d({0,0,0}) };
		 for(auto i=0;i<4;i++)
			 for(auto d=0; d<3;d++)
            J_[t][d] += refs_grad_[t][i] * var_.value(12 * t + 3 * i + d);

        det_[t] = J_[t].det();
        detmin_ = std::min(detmin_, det_[t]);
        ninverted_ += (det_[t] <= 0);

        K_[t] = dual_basis(J_[t]);
    }
    if (options_.debug > 3) std::cerr << "ok" << std::endl;

}
/*----------------------------------------------------------------------------*/
double Elliptic_smoother_3D::evaluate_energy() {
    evaluate_jacobian();
    double E = 0;
	 for(auto t=0; t<N_tets_;t++) {
        double c = chi(eps_, det_[t]);
        double f = (J_[t][0].dot(J_[t][0])+ J_[t][1].dot(J_[t][1]) + J_[t][2].dot(J_[t][2])) / std::pow(c, 2. / 3.);
        double g = (1 + det_[t] * det_[t]) / c;
        E += ((1 - options_.theta) * f + options_.theta * g) * vol_[t];
    }
    return E;
}
/*----------------------------------------------------------------------------*/
void Elliptic_smoother_3D::standard_elliptic_energy(const std::vector<double>&, double& F, std::vector<double>& G) {
    F += evaluate_energy();
	 for(auto t=0; t<N_tets_;t++) {
        double c1 = chi(eps_, det_[t]);
        double c2 = std::pow(c1, 2. / 3.);
        double c3 = chi_deriv(eps_, det_[t]);
        double f = (J_[t][0].dot(J_[t][0]) + J_[t][1].dot(J_[t][1]) + J_[t][2].dot(J_[t][2])) / c2;
        double g = (1 + det_[t] * det_[t]) / c1;

		  for(auto d=0; d<3; d++) {
			  math::Vector3d dfda = J_[t][d] * (2. / c2) - K_[t][d] * ((2. * f * c3) / (3. * c1));
			  math::Vector3d dgda = K_[t][d] * ((2 * det_[t] - g * c3) / c1);

			  for(auto tc=0; tc<4; tc++) {
				  double gradi = (dfda * (1. - options_.theta) + dgda * options_.theta).dot(refs_grad_[t][tc]) * vol_[t];
				  var_.add2grad(12 * t + 3 * tc + d, gradi, G);
			  }
		  }
	 }
}
/*----------------------------------------------------------------------------*/
void Elliptic_smoother_3D::standard_elliptic_energy_w_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<EllipticSmoothingMeshGlue::sparse_term>& H) {
    F += evaluate_energy();
	 for(auto t=0; t<N_tets_;t++) {
        double c1 = chi(eps_, det_[t]);
        double c2 = std::pow(c1, 2. / 3.);
        double c3 = chi_deriv(eps_, det_[t]);
        double f = (J_[t][0].dot(J_[t][0]) + J_[t][1].dot(J_[t][1]) + J_[t][2].dot(J_[t][2])) / c2;
        double g = (1 + det_[t] * det_[t]) / c1;


		  for(auto d=0; d<3;d++) {
			  math::Vector3d dfda = J_[t][d] * (2. / c2) - K_[t][d] * ((2. * f * c3) / (3. * c1));
			  math::Vector3d dgda = K_[t][d] * ((2 * det_[t] - g * c3) / c1);
            for(auto tc=0; tc<4; tc++) {
                double gradi = (dfda * (1. - options_.theta) + dgda * options_.theta) .dot(refs_grad_[t][tc]) * vol_[t];
                var_.add2grad(12 * t + 3 * tc + d, gradi, G);
            }
        }
    }

}
/*----------------------------------------------------------------------------*/
bool Elliptic_smoother_3D::run_lbfgs(std::vector<double>& X) {
  /*  const STLBFGS::Optimizer::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
        std::fill(G.begin(), G.end(), 0);
        F = 0;

        var_.set_reducted_values(X);
        energy(X, F, G);
    };
    STLBFGS::Optimizer opt = { func };
*/
 /*   opt.gtol = options_.bfgs_threshold;
    opt.ftol = options_.bfgs_threshold;
    opt.maxiter = options_.bfgs_maxiter;
    opt.verbose = options_.debug;
    opt.invH.history_depth = 5;*/
    //UM::LBFGS_Optimizer opt(func);
    //opt.gtol = options_.bfgs_threshold;
    //opt.maxiter = options_.bfgs_maxiter;
    //opt.verbose = options_.debug;
  //  opt.run(X);
    return true;
}
/*----------------------------------------------------------------------------*/
bool Elliptic_smoother_3D::run_newton(std::vector<double>& X) {
    return true;
}
/*----------------------------------------------------------------------------*/
bool Elliptic_smoother_3D::go() {
    if (options_.debug > 0) {
        std::cerr << "==== Running Elliptic smoother 3d. ====" << "\n";
        display_options(options_);
    }
    if (var_.get_nb_reduced_variables() == 0) {
        std::cerr << "No variables to optimize" << "\n";
        return true;
    }
    double total_vol_ = 0; // so that energy is always roughly the same
	 for(auto i=0; i<N_tets_; i++)
		 total_vol_ += vol_[i];
    if (options_.debug > 0) std::cerr << "total vol = " << total_vol_ << "\n";
	 for(auto i=0; i<N_tets_; i++)
		 vol_[i] = vol_[i] / total_vol_;
    evaluate_jacobian();
    double e0 = 1e-3;
    if (options_.eps_from_theorem) {
        eps_ = start_eps;
        double sigma = 0.5;
        double mu = (1 - sigma) * chi(eps_, detmin_);
        if (detmin_ < mu)
            eps_ = 2 * std::sqrt(mu * (mu - detmin_));
        else eps_ = 1e-10;
    }
    if (options_.barrier) {
        e0 = 1e-7;
        eps_ = 1e-7;
    }
    for(auto iter=0; iter<options_.maxiter; iter++) {
        if (options_.debug > 0) std::cerr << "iteration #" << iter << "\n";
        if (!options_.eps_from_theorem) {
            if (iter && iter % 10 == 0 && e0 > 1e-10) e0 /= 2.;
            eps_ = detmin_ > 0 ? e0 : std::sqrt(e0 * e0 + 0.04 * detmin_ * detmin_);
        }
        if (options_.debug > 0) std::cerr << "E: " << evaluate_energy() << " eps: " << eps_ << " detmin: " << detmin_ << " ninv: " << ninverted_ << std::endl;


        std::vector<double> X(var_.get_nb_reduced_variables());
		  var_.get_reduced_values(X);
        double E_prev = evaluate_energy();
        if (options_.use_newton) run_newton(X);
        else run_lbfgs(X);
		  var_.set_reduced_values(X);

        double E = evaluate_energy();
        if (options_.debug > 0) std::cerr << " E : " << E_prev << "  --->  " << E << "\n";
        if (options_.eps_from_theorem) {
            double sigma = std::max(1. - E / E_prev, 0.1);
            double mu = (1 - sigma) * chi(eps_, detmin_);
            if (detmin_ < mu)
                eps_ = 2 * std::sqrt(mu * (mu - detmin_));
            else eps_ = 1e-10;
        }
        iter_call_back(iter);
        if (detmin_ > 0 && std::abs(E_prev - E) / E < options_.static_threshold) break;
        if (options_.stopping_when_static && std::abs(E_prev - E) / E < options_.static_threshold) {
            std::cerr << "Stopping as Energy is static with threshold: " << options_.static_threshold << std::endl;
            break;
        }

    }
    if (options_.debug > 0) std::cerr << "E: " << evaluate_energy() << " detmin: " << detmin_ << " ninv: " << ninverted_ << std::endl;
    return !ninverted_;
}
/*----------------------------------------------------------------------------*/
