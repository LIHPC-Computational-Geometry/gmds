/*----------------------------------------------------------------------------*/
#  ifndef GMDS_ELLIPTIC_SMOOTHING_H
#	define GMDS_ELLIPTIC_SMOOTHING_H
/*----------------------------------------------------------------------------*/
#	include <gmds/math/Matrix.h>
#	include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
#	include <array>
#	include <assert.h>
#	include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds {
namespace smoothy {
/*----------------------------------------------------------------------------*/
class EllipticSmoothingMeshGlue
{
 public:
	/** return ne number of reduced variables */
	virtual int get_nb_reduced_variables() = 0;
	/** set the reduced values
	 * */
	virtual void set_reduced_values(const std::vector<double> &AX) = 0;
	/** return ne number of reduced variables */
	virtual void get_reduced_values(std::vector<double> &AX) = 0;

	virtual double value(const int i) = 0;
	virtual void add2grad(const int i, const double value, std::vector<double> &G) = 0;
	// Do not need to be unique in H
	typedef std::pair<std::array<int, 2>, double> sparse_term;
	virtual void add2Hessian(const int i, const int j, const double value, std::vector<sparse_term> &H) = 0;
};
/*----------------------------------------------------------------------------*/
struct EllipticSmoothingOptions
{
	EllipticSmoothingOptions(double _theta,
	                         int _maxiter,
	                         double _bfgs_threshold,
	                         int _bfgs_maxiter,
	                         int _debug,
	                         bool eps_from_theorem_,
	                         bool stopping_when_static_,
	                         double static_threshold_,
	                         bool _barrier,
	                         bool _use_newton) :
	  theta(_theta),
	  maxiter(_maxiter),
	  bfgs_threshold(_bfgs_threshold),
	  bfgs_maxiter(_bfgs_maxiter),
	  debug(_debug),
	  eps_from_theorem(eps_from_theorem_),
	  stopping_when_static(stopping_when_static_),
	  static_threshold(static_threshold_),
	  barrier(_barrier),
	  use_newton(_use_newton)
	{
	}
	double theta = 1. / 128;
	int maxiter = 10000;
	double bfgs_threshold = .1;
	int bfgs_maxiter = 30000;
	int debug = 1;
	bool eps_from_theorem = false;
	bool stopping_when_static = false;
	double static_threshold = 1e-5;
	bool barrier = false;
	bool use_newton = false;
};
/*----------------------------------------------------------------------------*/
/**
 * Default smoothing options
 */
const EllipticSmoothingOptions ESO_default2D(1. / 2, 10000, 1e-4, 30000, 1, false, false, 1e-5, false, false);
const EllipticSmoothingOptions ESO_default3D(1. / 2., 10000, 1e-4, 300, 1, false, false, 1e-5, false, false);
/*----------------------------------------------------------------------------*/
class EllipticSmoothing2D
{
 public:
	EllipticSmoothing2D(EllipticSmoothingMeshGlue &var, const int nb_tri, const EllipticSmoothingOptions &options = ESO_default2D) :
	  var_(var), m_nb_tri(nb_tri), m_options(options)
	{
		init();
	}
	// v1x, v1y, v2x, v2y, v3x, v3y
	void set_triangle_ref(const int t, const std::array<double, 6> &ref);
	void set_triangle_ref(const std::vector<std::array<double, 6>> &refs);
	void set_area(const int t, const double area);
	bool execute();

 public:
	virtual inline void energy(const std::vector<double> &X, double &F, std::vector<double> &G)
	{
		standard_elliptic_energy(X, F, G);
	}

	virtual inline void compute_energy_hessian(const std::vector<double> &X,
	                                           double &F,
	                                           std::vector<double> &G,
	                                           std::vector<EllipticSmoothingMeshGlue::sparse_term> &H)
	{
		standard_elliptic_energy_w_hessian(X, F, G, H);
	}

	void standard_elliptic_energy(const std::vector<double> &X, double &F, std::vector<double> &G);
	void standard_elliptic_energy_w_hessian(const std::vector<double> &X,
	                                        double &F,
	                                        std::vector<double> &G,
	                                        std::vector<EllipticSmoothingMeshGlue::sparse_term> &H);
	void init();

	typedef std::array<double, 6> tri_grad;
	typedef std::array<double, 4> mat22;
	virtual double evaluate_energy();
	void evaluate_jacobian();

	bool run_lbfgs(std::vector<double> &X);
	bool run_newton(std::vector<double> &X);

 private:
	// inputs
	EllipticSmoothingMeshGlue &var_;
	int m_nb_tri;
	EllipticSmoothingOptions m_options;
	std::vector<tri_grad> refs_grad_;

	// intern variables
	std::vector<mat22> m_J;         // per-tet Jacobian matrix = [[JX.x JX.y, JX.z], [JY.x, JY.y, JY.z], [JZ.x, JZ.y, JZ.z]]
	std::vector<mat22> m_K;         // per-tet dual basis: det J = dot J[i] * K[i]
	std::vector<double> m_det;      // per-tet determinant of the Jacobian matrix
	std::vector<double> m_area;     // per-tet determinant of the Jacobian matrix
	double m_eps;                   // regularization parameter, depends on min(jacobian)

	double m_detmin;       // min(jacobian) over all tetrahedra
	int m_nb_inverted;     // number of inverted tetrahedra

	double start_eps = 1.;
};
/*----------------------------------------------------------------------------*/
class Elliptic_smoother_3D
{
 public:
	Elliptic_smoother_3D(EllipticSmoothingMeshGlue &var, const int nb_tets, const EllipticSmoothingOptions &options = ESO_default3D) :
	  var_(var), N_tets_(nb_tets), options_(options)
	{
		init();
	}

#	ifdef WITH_UM
	void set_tet_ref(const int t, const std::array<UM::vec3, 4> tet_ref);
	void set_tet_ref(const std::vector<std::array<UM::vec3, 4>> tets_ref);
#	else
	void set_tet_ref(const int t, const std::array<math::Vector3d, 4> tet_ref);
	void set_tet_ref(const std::vector<std::array<math::Vector3d, 4>> tets_ref);
#	endif

	bool go();

 public:
	virtual inline void energy(const std::vector<double> &X, double &F, std::vector<double> &G)
	{
		standard_elliptic_energy(X, F, G);
	}
	virtual inline void iter_call_back(int iter_nb)
	{
		(void) iter_nb;
	}
	virtual inline void compute_energy_hessian(const std::vector<double> &X,
	                                           double &F,
	                                           std::vector<double> &G,
	                                           std::vector<EllipticSmoothingMeshGlue::sparse_term> &H)
	{
		standard_elliptic_energy_w_hessian(X, F, G, H);
	}

	void standard_elliptic_energy(const std::vector<double> &X, double &F, std::vector<double> &G);
	void standard_elliptic_energy_w_hessian(const std::vector<double> &X,
	                                        double &F,
	                                        std::vector<double> &G,
	                                        std::vector<EllipticSmoothingMeshGlue::sparse_term> &H);
	void init();

	double evaluate_energy();
	void evaluate_jacobian();

	bool run_lbfgs(std::vector<double> &X);
	bool run_newton(std::vector<double> &X);

	// inputs
	EllipticSmoothingMeshGlue &var_;
	int N_tets_;
	const EllipticSmoothingOptions options_;

#	ifdef WITH_UM
	std::vector<std::array<UM::vec3, 4>> refs_grad_;
#	else
	std::vector<std::array<math::Vector3d, 4>> refs_grad_;
#	endif

	// intern variables
	std::vector<math::Matrix33> J_;     // per-tet Jacobian matrix = [[JX.x JX.y, JX.z], [JY.x, JY.y, JY.z], [JZ.x, JZ.y, JZ.z]]
	std::vector<math::Matrix33> K_;     // per-tet dual basis: det J = dot J[i] * K[i]
	std::vector<double> det_;           // per-tet determinant of the Jacobian matrix
	std::vector<double> vol_;           // per-tet determinant of the Jacobian matrix
	double eps_;                        // regularization parameter, depends on min(jacobian)

	double detmin_;     // min(jacobian) over all tetrahedra
	int ninverted_;     // number of inverted tetrahedra

	double start_eps = 1.;
};
/*----------------------------------------------------------------------------*/
class EllipticSmoothingMeshGlue2D : public EllipticSmoothingMeshGlue
{
 public:
	EllipticSmoothingMeshGlue2D(const std::vector<double> &verts, const std::vector<std::array<int, 3>> &triangles, const std::vector<bool> &locks)
	{
		m_nb_unlocked = 0;
		m_vertices.resize(verts.size());
		m_locked.resize(verts.size());

		for (int i = 0; i < (int) verts.size(); i++)
			m_vertices[i] = verts[i];
		for (int i = 0; i < (int) verts.size(); i++)
			m_locked[i] = locks[i];

		m_map_to_vert.resize(6 * triangles.size());
		m_reduced.resize(m_vertices.size());
		for (int t = 0; t < (int) triangles.size(); t++) {
			for (int tv = 0; tv < 3; tv++) {
				m_map_to_vert[6 * t + 2 * tv + 0] = 2 * triangles[t][tv] + 0;
				m_map_to_vert[6 * t + 2 * tv + 1] = 2 * triangles[t][tv] + 1;
			}
		}
		for (int i = 0; i < (int) m_vertices.size(); i++) {
			if (!m_locked[i]) {
				m_reduced[i] = m_nb_unlocked++;
			}
			else
				m_reduced[i] = -1;
		}
	}
	/*----------------------------------------------------------------------------*/
	inline int get_nb_reduced_variables()
	{
		return m_nb_unlocked;
	}
	/*----------------------------------------------------------------------------*/
	inline void set_reduced_values(const std::vector<double> &X)
	{
		for (int i = 0; i < (int) m_vertices.size(); i++)
			if (!m_locked[i]) {
				m_vertices[i] = X[m_reduced[i]];
			}
	}
	inline void get_reduced_values(std::vector<double> &X)
	{
		for (int i = 0; i < (int) m_vertices.size(); i++)
			if (!m_locked[i]) {
				X[m_reduced[i]] = m_vertices[i];
			}
	}

	inline double value(const int i)
	{
		return m_vertices[m_map_to_vert[i]];
	}
	inline void add2grad(const int i, const double value, std::vector<double> &G)
	{
		if (!m_locked[m_map_to_vert[i]]) G[m_reduced[m_map_to_vert[i]]] += value;
	}
	inline void add2Hessian(const int i, const int j, const double value, std::vector<sparse_term> &H)
	{
		H.push_back({{m_reduced[m_map_to_vert[i]], m_reduced[m_map_to_vert[j]]}, value});
	}
	void get_verts(std::vector<double> &AVertices)
	{
		AVertices.resize(m_vertices.size());
		for (int i = 0; i < (int) m_vertices.size(); i++) {
			AVertices[i] = m_vertices[i];
		}
	}

 private:
	int m_nb_unlocked;
	std::vector<int> m_map_to_vert;
	std::vector<double> m_vertices;
	std::vector<bool> m_locked;
	std::vector<int> m_reduced;
};
/*----------------------------------------------------------------------------*/
class EllipticSmoothingMeshGlue3D : public EllipticSmoothingMeshGlue
{
 public:
	EllipticSmoothingMeshGlue3D(const std::vector<double> &verts, const std::vector<std::array<int, 4>> &tets, const std::vector<bool> &locks)
	{
		m_nb_unlocked = 0;
		m_verts.resize(verts.size());
		m_locked.resize(verts.size());

		for (int i = 0; i < (int) verts.size(); i++)
			m_verts[i] = verts[i];
		for (int i = 0; i < (int) verts.size(); i++)
			m_locked[i] = locks[i];

		m_map_to_vert.resize(12 * tets.size());
		m_reduc.resize(m_verts.size());
		for (int t = 0; t < (int) tets.size(); t++) {
			for (int tv = 0; tv < 4; tv++) {
				for (int d = 0; d < 3; d++) {
					m_map_to_vert[12 * t + 3 * tv + d] = 3 * tets[t][tv] + d;
				}
			}
		}
		for (int i = 0; i < (int) m_verts.size(); i++) {
			if (!m_locked[i]) {
				m_reduc[i] = m_nb_unlocked++;
			}
			else
				m_reduc[i] = -1;
		}
	}
	inline int get_nb_reduced_variables()
	{
		return m_nb_unlocked;
	}
	inline void set_reduced_values(const std::vector<double> &X)
	{
		for (int i = 0; i < (int) m_verts.size(); i++)
			if (!m_locked[i]) {
				m_verts[i] = X[m_reduc[i]];
			}
	}
	inline void get_reduced_values(std::vector<double> &X)
	{
		for (int i = 0; i < (int) m_verts.size(); i++)
			if (!m_locked[i]) {
				X[m_reduc[i]] = m_verts[i];
			}
	}
	inline double value(const int i)
	{
		return m_verts[m_map_to_vert[i]];
	}
	inline void add2grad(const int i, const double value, std::vector<double> &G)
	{
		if (!m_locked[m_map_to_vert[i]]) G[m_reduc[m_map_to_vert[i]]] += value;
	}
	inline void add2Hessian(const int i, const int j, const double value, std::vector<sparse_term> &H)
	{
		H.push_back({{m_reduc[m_map_to_vert[i]], m_reduc[m_map_to_vert[j]]}, value});
	}
	void get_verts(std::vector<double> &verts)
	{
		for (int i = 0; i < (int) m_verts.size(); i++) {
			verts[i] = m_verts[i];
		}
	}

 private:
	int m_nb_unlocked;
	std::vector<int> m_map_to_vert;
	std::vector<double> m_verts;
	std::vector<bool> m_locked;
	std::vector<int> m_reduc;
};
}     // namespace smoothy
}     // namespace gmds

/*----------------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------------*/
