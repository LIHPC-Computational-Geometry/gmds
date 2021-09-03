/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*---------------------------------------------------------------------------*/
// STL Header files
#include <iostream>
#include <map>
#include <set>
/*---------------------------------------------------------------------------*/
#include "FrameFieldSmoother.h"
/*---------------------------------------------------------------------------*/
// GMDS Header files
#include "GMDS/Math/Matrix.h"
/*---------------------------------------------------------------------------*/
// HLBFGS Header file
#include "HLBFGS.h"
#include <GMDS/IO/VTKWriter.h>

/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
// GLOBAL VARIABLES USED FOR THE CALL TO HLBFGS
/*---------------------------------------------------------------------------*/
gmds::IGMesh* HLBFGS_mesh = 0;
// mark to identify the nodes that are classified on points
int HLBFGS_mark_point = 0;
// mark to identify the nodes that are classified on curves
int HLBFGS_mark_curve = 0;
// mark to identify the nodes that are classified on surfaces
int HLBFGS_mark_surface = 0;
//the ordered list of candidate nodes
std::vector<gmds::Node> HLBFGS_candidates;
//the list of edges that connect candidates
// these edges will be used to compute smoothing values between
// end nodes
std::vector<gmds::Edge> HLBFGS_edges;
//the index of a candidate node in m_candidates
std::map<gmds::TCellID, int> HLBFGS_candidates_index;
//We keep in mind a Chart for each boundary node we work on
std::map<gmds::TCellID, gmds::math::Chart> HLBFGS_boundary_triad;

/*---------------------------------------------------------------------------*/
void initMatrix(Node& ANode,
		gmds::math::Matrix<3, 3,double>& Mix,
		gmds::math::Matrix<3, 3,double>& Miy,
		gmds::math::Matrix<3, 3,double>& Miz,
		gmds::math::Matrix<3, 3,double>& DMix,
		gmds::math::Matrix<3, 3,double>& DMiy,
		gmds::math::Matrix<3, 3,double>& DMiz,
		std::vector<TCoord>& trig_cos_vectors,
		std::vector<TCoord>& trig_sin_vectors)
{
  Node n = ANode;
  int index = HLBFGS_candidates_index[n.getID()];
  // contribution to 0-indexed vertex
  if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_curve))
    {
      //on a sharp edge
      Mix.set(0, 0, 1);
      Mix.set(1, 1, trig_cos_vectors[3 * index]);
      Mix.set(2, 2, trig_cos_vectors[3 * index]);
      Mix.set(1, 2, -trig_sin_vectors[3 * index]);
      Mix.set(2, 1, trig_sin_vectors[3 * index]);

      Miy.set(1, 1, 1);
      Miy.set(0, 0, trig_cos_vectors[3 * index + 1]);
      Miy.set(2, 2, trig_cos_vectors[3 * index + 1]);
      Miy.set(0, 2, trig_sin_vectors[3 * index + 1]);
      Miy.set(2, 0, -trig_sin_vectors[3 * index + 1]);

      Miz.set(2, 2, 1);
      Miz.set(0, 0, trig_cos_vectors[3 * index + 2]);
      Miz.set(1, 1, trig_cos_vectors[3 * index + 2]);
      Miz.set(0, 1, -trig_sin_vectors[3 * index + 2]);
      Miz.set(1, 0, trig_sin_vectors[3 * index + 2]);

      // Steady values for the computation, so zero derivatives
    }
  else if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_surface))
    {
      Mix.set(0, 0, 1.0);
      Mix.set(1, 1, 1.0);
      Mix.set(2, 2, 1.0);

      Miy.set(0, 0, 1.0);
      Miy.set(1, 1, 1.0);
      Miy.set(2, 2, 1.0);

      double VX[3];
      double VY[3];
      math::Chart t = HLBFGS_boundary_triad[n.getID()];
      math::Vector tX = t.X();
      math::Vector tY = t.Y();
      math::Vector tZ = t.Z();

      VX[0] = trig_cos_vectors[3 * index + 2] * tY[0] + trig_sin_vectors[3 * index + 2] * tZ[0];
      VX[1] = trig_cos_vectors[3 * index + 2] * tY[1] + trig_sin_vectors[3 * index + 2] * tZ[1];
      VX[2] = trig_cos_vectors[3 * index + 2] * tY[2] + trig_sin_vectors[3 * index + 2] * tZ[2];

      VY[0] = -trig_sin_vectors[3 * index + 2] * tY[0] + trig_cos_vectors[3 * index + 2] * tZ[0];
      VY[1] = -trig_sin_vectors[3 * index + 2] * tY[1] + trig_cos_vectors[3 * index + 2] * tZ[1];
      VY[2] = -trig_sin_vectors[3 * index + 2] * tY[2] + trig_cos_vectors[3 * index + 2] * tZ[2];

      Miz.set(0, 0, VX[0]);
      Miz.set(1, 0, VX[1]);
      Miz.set(2, 0, VX[2]);

      Miz.set(0, 1, VY[0]);
      Miz.set(1, 1, VY[1]);
      Miz.set(2, 1, VY[2]);

      Miz.set(0, 2, tX[0]);
      Miz.set(1, 2, tX[1]);
      Miz.set(2, 2, tX[2]);

      DMiz.set(0, 0, -trig_sin_vectors[3 * index + 2] * tY[0] + trig_cos_vectors[3 * index + 2] * tZ[0]);
      DMiz.set(1, 0, -trig_sin_vectors[3 * index + 2] * tY[1] + trig_cos_vectors[3 * index + 2] * tZ[1]);
      DMiz.set(2, 0, -trig_sin_vectors[3 * index + 2] * tY[2] + trig_cos_vectors[3 * index + 2] * tZ[2]);

      DMiz.set(0, 1, -trig_cos_vectors[3 * index + 2] * tY[0] - trig_sin_vectors[3 * index + 2] * tZ[0]);
      DMiz.set(1, 1, -trig_cos_vectors[3 * index + 2] * tY[1] - trig_sin_vectors[3 * index + 2] * tZ[1]);
      DMiz.set(2, 1, -trig_cos_vectors[3 * index + 2] * tY[2] - trig_sin_vectors[3 * index + 2] * tZ[2]);

      DMiz.set(0, 2, 0.0);
      DMiz.set(1, 2, 0.0);
      DMiz.set(2, 2, 0.0);

    }
  else{ //General case
    //We are inside the domain and not on the boundary, so traditionnal matrix filling in
    //std::cout<<"vi inside"<<std::endl;
    Mix.set(0, 0, 1);
    Mix.set(1, 1, trig_cos_vectors[3 * index]);
    Mix.set(2, 2, trig_cos_vectors[3 * index]);
    Mix.set(1, 2, -trig_sin_vectors[3 * index]);
    Mix.set(2, 1, trig_sin_vectors[3 * index]);

    DMix.set(1, 1, -trig_sin_vectors[3 * index]);
    DMix.set(2, 2, -trig_sin_vectors[3 * index]);
    DMix.set(1, 2, -trig_cos_vectors[3 * index]);
    DMix.set(2, 1, trig_cos_vectors[3 * index]);

    Miy.set(1, 1, 1);
    Miy.set(0, 0, trig_cos_vectors[3 * index + 1]);
    Miy.set(2, 2, trig_cos_vectors[3 * index + 1]);
    Miy.set(0, 2, trig_sin_vectors[3 * index + 1]);
    Miy.set(2, 0, -trig_sin_vectors[3 * index + 1]);

    DMiy.set(0, 0, -trig_sin_vectors[3 * index + 1]);
    DMiy.set(2, 2, -trig_sin_vectors[3 * index + 1]);
    DMiy.set(0, 2, trig_cos_vectors[3 * index + 1]);
    DMiy.set(2, 0, -trig_cos_vectors[3 * index + 1]);

    Miz.set(2, 2, 1);
    Miz.set(0, 0, trig_cos_vectors[3 * index + 2]);
    Miz.set(1, 1, trig_cos_vectors[3 * index + 2]);
    Miz.set(0, 1, -trig_sin_vectors[3 * index + 2]);
    Miz.set(1, 0, trig_sin_vectors[3 * index + 2]);

    DMiz.set(0, 0, -trig_sin_vectors[3 * index + 2]);
    DMiz.set(1, 1, -trig_sin_vectors[3 * index + 2]);
    DMiz.set(0, 1, -trig_cos_vectors[3 * index + 2]);
    DMiz.set(1, 0, trig_cos_vectors[3 * index + 2]);
  }
}
/*----------------------------------------------------------------------------*/
/*
 * IN N -> nb variables to consider (3xnb vertices in our case)
 * IN x -> size N contient les angles d'Euler pour chaque sommet à traiter
 *  		initialisé avant
 * IN prev_x -> NON UTILISE valeur précédente de x sinon
 * OUT func-> pointer on the computed function
 * OUT grad -> grad(func)
 */
void evalF(int N, double* x, double *prev_x, double* func, double* grad)
{
  size_t nvars = N;

  //Size of "grad" is already equal to N, since it points on a vector
  //initialized with the right size in HLBFGS

  //Initialization of trigo function for incoming differetial computations
  std::vector<TCoord> trig_cos_vectors, trig_sin_vectors;
  trig_cos_vectors.resize(N);
  trig_sin_vectors.resize(N);

  for (unsigned int i = 0; i < N; i++)
    {
      //std::cout << "X[" << i << "] = " << x[i] << std::endl;
      trig_cos_vectors[i] = cos(x[i]);
      trig_sin_vectors[i] = sin(x[i]);
    }

  *func = 0.0;

  // We do not initialze grad since nothing said
  for (unsigned int i = 0; i < N; i++)
    grad[i] = 0.0;


  // the matrix A to solve is build by going throug all the mesh edges
  for (unsigned int i_edge = 0; i_edge < HLBFGS_edges.size(); i_edge++) {
    std::vector<Node> edge_nodes = HLBFGS_edges[i_edge].get<Node>();
    Node ni = edge_nodes[0];
    Node nj = edge_nodes[1];
      
    int indexi = HLBFGS_candidates_index[ni.getID()];
    int indexj = HLBFGS_candidates_index[nj.getID()];


    gmds::math::Matrix<3, 3,double> Mix, Miy, Miz;
    gmds::math::Matrix<3, 3,double> Mjx, Mjy, Mjz;
    gmds::math::Matrix<3, 3,double> DMix, DMiy, DMiz;
    gmds::math::Matrix<3, 3,double> DMjx, DMjy, DMjz;
    gmds::math::Matrix<3, 3,double> M, Mi, Mj, Mi1, Mi2, Mi3;
    // Matrix initialization. It is mandatory to avoid numerical issues
    // depending on the used compiler
    TCoord diff1[3][3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    TCoord diff2[3][3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    TCoord diffi[3] = { 0, 0, 0 };
    TCoord diffj[3] = { 0, 0, 0 };
    TCoord e[18];

    // Matrix computation for node n0
    initMatrix(ni, Mix, Miy, Miz, DMix, DMiy, DMiz, trig_cos_vectors, trig_sin_vectors);

    Mi = Mix * Miy * Miz;
    Mi1 = DMix * Miy * Miz;
    Mi2 = Mix * DMiy * Miz;
    Mi3 = Mix * Miy * DMiz;

    // Matrix computation for node n1
    initMatrix(nj, Mjx, Mjy, Mjz, DMjx, DMjy, DMjz, trig_cos_vectors, trig_sin_vectors);

    Mj = Mjx * Mjy * Mjz;

    //WARNING, I could make a mistake and transpose the wrong Matrix
    //(idem for Mi1, Mi2, Mi3)
    M = Mi.transpose() * Mj;
    //M = Mi * Mj.transpose();
    //M = Mj.transpose() * Mi;

    gmds::math::Matrix<3, 3, double> tmpMat = (Mi.transpose() * (DMjx * Mjy * Mjz));

    tmpMat.getTab(diff2[0]);

    tmpMat = (Mi.transpose() * (Mjx * DMjy * Mjz));
    tmpMat.getTab(diff2[1]);


    tmpMat = (Mi.transpose() * (Mjx * Mjy * DMjz));
    tmpMat.getTab(diff2[2]);

    (Mi1.transpose() * Mj).getTab(diff1[0]);
    (Mi2.transpose() * Mj).getTab(diff1[1]);
    (Mi3.transpose() * Mj).getTab(diff1[2]);
    //OR 
    //		(Mj.transpose() * Mi1).getTab(diff1[0]);
    //		(Mj.transpose() * Mi2).getTab(diff1[1]);
    //		(Mj.transpose() * Mi3).getTab(diff1[2]);
    //OR
    //		(Mi1 * Mj.transpose()).getTab(diff1[0]);
    //		(Mi2 * Mj.transpose()).getTab(diff1[1]);
    //		(Mi3 * Mj.transpose()).getTab(diff1[2]);


    e[0] = M.get(0, 0) * M.get(0, 1); 
    e[1] = M.get(0, 0) * M.get(0, 2); 
    e[2] = M.get(0, 1) * M.get(0, 2);
    e[3] = M.get(1, 0) * M.get(1, 1); 
    e[4] = M.get(1, 0) * M.get(1, 2); 
    e[5] = M.get(1, 1) * M.get(1, 2);
    e[6] = M.get(2, 0) * M.get(2, 1); 
    e[7] = M.get(2, 0) * M.get(2, 2); 
    e[8] = M.get(2, 1) * M.get(2, 2);
    e[9] = M.get(0, 0) * M.get(1, 0); 
    e[10] = M.get(0, 0) * M.get(2, 0); 
    e[11] = M.get(1, 0) * M.get(2, 0);
    e[12] = M.get(0, 1) * M.get(1, 1); 
    e[13] = M.get(0, 1) * M.get(2, 1); 
    e[14] = M.get(1, 1) * M.get(2, 1);
    e[15] = M.get(0, 2) * M.get(1, 2); 
    e[16] = M.get(0, 2) * M.get(2, 2); 
    e[17] = M.get(1, 2) * M.get(2, 2);

    double sigma = 1.0;

    double local_f = 0.0;


    static const int eflag[9][2][2] = {
      { { 0, 0 }, { 0, 1 } }, { { 0, 0 }, { 0, 2 } }, { { 0, 1 }, { 0, 2 } },
      { { 1, 0 }, { 1, 1 } }, { { 1, 0 }, { 1, 2 } }, { { 1, 1 }, { 1, 2 } },
      { { 2, 0 }, { 2, 1 } }, { { 2, 0 }, { 2, 2 } }, { { 2, 1 }, { 2, 2 } } };


    for (int k = 0; k < 9; k++)
      {
	double v2 = exp(e[k + 9] * e[k + 9] / sigma);

	local_f += 2 * (v2 - 1);
	for (int l = 0; l < 3; l++)
	  {
	    diffi[l] += 2 * (e[k + 9] * (diff1[l][eflag[k][0][1]][eflag[k][0][0]] * M.get(eflag[k][1][1], eflag[k][1][0]) + 
					 diff1[l][eflag[k][1][1]][eflag[k][1][0]] * M.get(eflag[k][0][1], eflag[k][0][0])) * v2 / sigma);
	    diffj[l] += 2 * (e[k + 9] * (diff2[l][eflag[k][0][1]][eflag[k][0][0]] * M.get(eflag[k][1][1], eflag[k][1][0]) + diff2[l][eflag[k][1][1]][eflag[k][1][0]] * M.get(eflag[k][0][1], eflag[k][0][0])) * v2 / sigma);
	  }
      }

    
    *func += (local_f / 2.0);
    for (int k = 0; k < 3; k++)
      {
	grad[3 * indexi + k] += diffi[k];
	grad[3 * indexj + k] += diffj[k];
      }

  }//  for (unsigned int i_edge = 0; i_edge < HLBFGS_edges.size(); i_edge++)


}
/*----------------------------------------------------------------------------*/
void newiteration(int iter, int call_iter, double *x, double* f, double *g, 
		  double* gnorm)
{
  std::cout << iter << ": " << call_iter << " " << *f << " " << *gnorm << std::endl;
}

/*----------------------------------------------------------------------------*/
void FrameFieldSmoother::rebuildQuaternions(double*& eulerAngles)
{

  for (unsigned int i = 0; i < HLBFGS_candidates.size(); i++)
    {


      Node n = HLBFGS_candidates[i];
      TCoord angleX = eulerAngles[3 * i];
      TCoord angleY = eulerAngles[3 * i + 1];
      TCoord angleZ = eulerAngles[3 * i + 2];

      math::Quaternion q;
      if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_curve))
	{
	  // on a curve
	  q.setFromEulerAngle(angleX, angleY, angleZ);
	}
      else if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_surface))
	{
	  // on a surface
			

	  math::Chart t = HLBFGS_boundary_triad[n.getID()];

	  math::Vector newVY = cos(angleZ) * t.Y();
	  math::Vector newVZ = sin(angleZ) * t.Z();
	  math::Vector tmp1 = newVY + newVZ;

	  math::Vector tmp2 = t.X().cross(tmp1);

	  q = math::Quaternion(math::Chart(t.X(), tmp1, tmp2));
	  //TODO : change if we put Z = normal
	}
      else {
	//inside the volume
	q.setFromEulerAngle(angleX, angleY, angleZ);
      }

      //Quaternion update
      if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_surface))
	{
	  math::Vector normal = (*m_surf_normal)[n.getID()];
	  q = q.alignWith(normal);
	}
      (*m_frame_field)[n.getID()] = q;

    }
}

/*---------------------------------------------------------------------------*/
FrameFieldSmoother::
FrameFieldSmoother(gmds::IGMesh* AMesh,
		   gmds::Variable<gmds::math::Quaternion>*& AField,
		   gmds::Variable<gmds::math::Vector>*& ANormal)
  : m_frame_field(AField), m_surf_normal(ANormal),
    m_mark_candidates(0), m_candidate_mark_initialized(false),
    m_boundary_marks_initialized(false)
{
  m_mesh = AMesh;
  HLBFGS_mesh = m_mesh;
}
/*---------------------------------------------------------------------------*/
void FrameFieldSmoother::selectNodes(const int AMark)
{
  m_mark_candidates = AMark;
  m_candidate_mark_initialized = true;
}
/*---------------------------------------------------------------------------*/
void FrameFieldSmoother::
initBoundaryMarks(const int AMarkPnt, const int AMarkCurve,
		  const int AMarkSurf)
{
  HLBFGS_mark_point = AMarkPnt;
  HLBFGS_mark_curve = AMarkCurve;
  HLBFGS_mark_surface = AMarkSurf;
  m_boundary_marks_initialized = true;
}
/*---------------------------------------------------------------------------*/
void FrameFieldSmoother::execute()
{
  if (!m_candidate_mark_initialized)
    throw GMDSException("Candidate mark is not initialized");
  if (!m_boundary_marks_initialized)
    throw GMDSException("Boundary marks are not initialized");

  performSmoothing();
}
/*---------------------------------------------------------------------------*/
void FrameFieldSmoother::initCandidates()
{
  HLBFGS_candidates.clear();
  IGMesh::node_iterator it_nodes = HLBFGS_mesh->nodes_begin();
  for (; !it_nodes.isDone(); it_nodes.next())
    {
      Node n = it_nodes.value();
      if (HLBFGS_mesh->isMarked(n, m_mark_candidates))
	HLBFGS_candidates.push_back(it_nodes.value());
    }
  std::cout << std::endl << "Nb candidates (" << HLBFGS_candidates.size()
	    << " / " << HLBFGS_mesh->getNbNodes() << ")" << std::endl;;
}
/*---------------------------------------------------------------------------*/
int FrameFieldSmoother::initOptimizationData(double*& eulerAngles)
{
  int nb_nodes_on_curves = 0;
  int nb_nodes_on_surfaces = 0;
  int nb_nodes_in_volumes = 0;

  // Euler angles are computed for every defined vertex (including ridges)
  for (unsigned int i = 0; i < HLBFGS_candidates.size(); i++)
    {
      Node n = HLBFGS_candidates[i];
      TCellID n_id = n.getID();
      //we store the candidate index of n_i
      HLBFGS_candidates_index[n_id] = i;

      math::Quaternion q = (*m_frame_field)[n_id];

      if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_curve))
	{ //node on an edge??
	  q.toEulerAngle(eulerAngles[3 * i],
			 eulerAngles[3 * i + 1],
			 eulerAngles[3 * i + 2]);
	  nb_nodes_on_curves++;
	}
      else if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_surface))
	{

	  nb_nodes_on_surfaces++;
	  math::Chart t(q);
	  math::Vector X[6];
	  X[0] = t.X();
	  X[1] = t.Y();
	  X[2] = t.Z();
	  //inv of X[0]
	  X[3][0] = -X[0][0];
	  X[3][1] = -X[0][1];
	  X[3][2] = -X[0][2];
	  //inv of X[1]
	  X[4][0] = -X[1][0];
	  X[4][1] = -X[1][1];
	  X[4][2] = -X[1][2];
	  //inv of X[2]
	  X[5][0] = -X[2][0];
	  X[5][1] = -X[2][1];
	  X[5][2] = -X[2][2];


	  math::Vector normal;
	  math::Vector normal_vec = (*m_surf_normal)[n_id];
	  normal[0] = normal_vec.X();
	  normal[1] = normal_vec.Y();
	  normal[2] = normal_vec.Z();

	  // We get among X, the vector that is the most aligned with
	  // the vertex normal
	  int indexNormal = 0;
	  TCoord maxVal = normal.dot(X[0]);
	  for (unsigned int j = 1; j<6; j++)
	    {
	      TCoord currentVal = normal.dot(X[j]);
	      if (currentVal>maxVal)
		{
		  indexNormal = j;
		  maxVal = currentVal;
		}
	    }

	  // On construit la Chart externe
	  math::Vector ChartRes[3];
	  ChartRes[0] = X[indexNormal];
	  if (indexNormal == 0){
	    ChartRes[1] = X[1]; ChartRes[2] = X[2];
	  }
	  else if (indexNormal == 1){
	    ChartRes[1] = X[2]; ChartRes[2] = X[0];
	  }
	  else if (indexNormal == 2){
	    ChartRes[1] = X[0]; ChartRes[2] = X[1];
	  }
	  else if (indexNormal == 3){
	    ChartRes[1] = X[2]; ChartRes[2] = X[1];
	  }
	  else if (indexNormal == 4){
	    ChartRes[1] = X[0]; ChartRes[2] = X[2];
	  }
	  else if (indexNormal == 5){
	    ChartRes[1] = X[0]; ChartRes[2] = X[4];
	  }

	  HLBFGS_boundary_triad[n_id] =
	    math::Chart(ChartRes[0], ChartRes[1], ChartRes[2]);
		
	  eulerAngles[3 * i] = 0;
	  eulerAngles[3 * i + 1] = 0;
	  eulerAngles[3 * i + 2] = 0;
        
    }
      else {//general case in the volume
          q.toEulerAngle(eulerAngles[3 * i],
                         eulerAngles[3 * i + 1],
                         eulerAngles[3 * i + 2]);
          
          nb_nodes_in_volumes++;
      }
    }
    
    //Definition of the support Edges
    HLBFGS_edges.clear();
    for (unsigned int i = 0; i < HLBFGS_candidates.size(); i++){
        Node n = HLBFGS_candidates[i];
        
        // this test must be useless, but it is to be sure (Debug purpose)
        if (!m_mesh->isMarked(n, m_mark_candidates))
            continue;
        
        std::vector<Edge> adj_edges = n.get<Edge>();
        
        for (unsigned int j = 0; j < adj_edges.size(); j++)
        {
            Edge e_j = adj_edges[j];
            
            std::vector<Node> e_j_nodes = e_j.get<Node>();
            Node other_node =
            (e_j_nodes[0].getID() == n.getID()) ? e_j_nodes[1] : e_j_nodes[0];
            
//            if (m_mesh->isMarked(other_node, m_mark_candidates))
            {
                if (other_node.getID() < n.getID()) //to store the edge only once
                    HLBFGS_edges.push_back(e_j);
            }
        }
  }

  std::cout << "Nb edges: " << HLBFGS_edges.size() << std::endl;
  return nb_nodes_in_volumes;
}
/*---------------------------------------------------------------------------*/
void FrameFieldSmoother::performSmoothing()
{
  //We collect the nodes we want to smooth the Quaternion on
  initCandidates();

  int euler_angles_size = 3 * HLBFGS_candidates.size();
  double* euler_angles = new double[euler_angles_size];

  int nb_candidates = initOptimizationData(euler_angles); 
  //new de euler_angle dans la fonction init
	


  std::cout << "Global Euler Smoothing with HLBFGS"
	    << " (nb candidates= " << HLBFGS_candidates.size() << ")" << std::endl;

  //Init. de HLBFGS
  double parameter[20];
  int info[20];
  int N = 3 * HLBFGS_candidates.size();// Nb variable, 3 euler angles per node
  int M = 7;//Default choice for the optimization algorithm
  int T = 0;
  int num_iter = 1000;
  bool with_hessian = false;
  INIT_HLBFGS(parameter, info);

  info[3] = 1; //0 is default value, but 1 is recommended is the HLBGS documentation
  info[4] = num_iter;
  info[6] = T;
  info[7] = with_hessian ? 1 : 0;
  info[10] = 0;
  info[11] = 1;


  std::cout << "=== Before HLBFGS Smoothing" << std::endl;

  int nannb = 0;
  for (unsigned int i = 0; i<N; i++){
    double d = euler_angles[i];
    if (!(d>-5) && !(d < 5)){
      nannb++;
    }

  }
  std::cout << "Nb NAN in Euler angles: " << nannb << std::endl;

  HLBFGS(N, M, &euler_angles[0], evalF, 0, HLBFGS_UPDATE_Hessian,
	 newiteration, parameter, info);
  
  std::cout << "=== After smoothing, nb iter : " << info[2]
	    << " for " << info[1] << "evaluations" <<std::endl;

  std::cout << "Quaternions rebuild" << std::endl;
  rebuildQuaternions(euler_angles);
  std::cout << "Quaternions rebuild done" << std::endl;


  delete[] euler_angles;
    
    
    
    
    static int nb_file = 0;
    
    IGMesh::node_iterator it = m_mesh->nodes_begin();
    double bound = 100000;
    double x_min = bound;
    double y_min = bound;
    double z_min = bound;
    double x_max = -bound;
    double y_max = -bound;
    double z_max = -bound;
    for (; !it.isDone(); it.next())
    {
        Node n = it.value();
        math::Point p = n.getPoint();
        if (p.X() < x_min)
            x_min = p.X();
        if (p.X() > x_max)
            x_max = p.X();
        
        if (p.Y() < y_min)
            y_min = p.Y();
        if (p.Y() > y_max)
            y_max = p.Y();
        
        if (p.Z() < z_min)
            z_min = p.Z();
        if (p.Z() > z_max)
            z_max = p.Z();
    }
    double dist_x = x_max - x_min;
    double dist_y = y_max - y_min;
    double dist_z = z_max - z_min;
    
    double cube_size = 0;
    if (dist_x <= dist_y && dist_x <= dist_z){
        cube_size = dist_x;
    }
    else if (dist_y <= dist_x && dist_y <= dist_z){
        cube_size = dist_y;
    }
    else
        cube_size = dist_z;
    
    //    cube_size /= 20;
    
    MeshModel model_cube(DIM3 | R | N | R2N);
    IGMesh mesh_cube(model_cube);
    
    Variable<int>* v = mesh_cube.newVariable<int>(GMDS_REGION, "classification");
    
    
    for (it = m_mesh->nodes_begin(); !it.isDone(); it.next())
    {
        Node n = it.value();
        math::Point center = n.getPoint();
        
        math::Quaternion q =(*m_frame_field)[n.getID()];
        math::Chart c(q);
        
        math::Vector3d evx(c.X()[0],c.X()[1],c.X()[2]);
        math::Vector3d evy(c.Y()[0],c.Y()[1],c.Y()[2]);
        math::Vector3d evz(c.Z()[0],c.Z()[1],c.Z()[2]);
        
        evx.normalize();
        evy.normalize();
        evz.normalize();
        
        if(evx.X()>100){
            std::cout<<"Error Nan for node "<<n.getID()<<std::endl;
            continue;
        }
        math::Vector vx(evx.X(),evx.Y(),evx.Z());
        math::Vector vy(evy.X(),evy.Y(),evy.Z());
        math::Vector vz(evz.X(),evz.Y(),evz.Z());
        
        math::Point p1 = center + (vx + vy - vz)*cube_size;
        
        Node n1 = mesh_cube.newNode(p1);
        math::Point p2 = center + (vx - vy - vz)*cube_size;
        Node n2 = mesh_cube.newNode(p2);
        math::Point p3 = center + (vx + vy + vz).opp()*cube_size;
        Node n3 = mesh_cube.newNode(p3);
        math::Point p4 = center + (vy - vx - vz)*cube_size;
        Node n4 = mesh_cube.newNode(p4);
        
        math::Point p5 = center + (vx + vy + vz)*cube_size;
        Node n5 = mesh_cube.newNode(p5);
        math::Point p6 = center + (vx - vy + vz)*cube_size;
        Node n6 = mesh_cube.newNode(p6);
        math::Point p7 = center + (vx + vy - vz).opp()*cube_size;
        Node n7 = mesh_cube.newNode(p7);
        math::Point p8 = center + (vy - vx + vz)*cube_size;
        Node n8 = mesh_cube.newNode(p8);
        //        std::cout<<"Cube from "<<std::endl;
        //        std::cout<<"\t "<<p1<<std::endl;
        //        std::cout<<"\t "<<p2<<std::endl;
        //        std::cout<<"\t "<<p3<<std::endl;
        //        std::cout<<"\t "<<p4<<std::endl;
        //        std::cout<<"\t "<<p5<<std::endl;
        //        std::cout<<"\t "<<p6<<std::endl;
        //        std::cout<<"\t "<<p7<<std::endl;
        //        std::cout<<"\t "<<p8<<std::endl;
        Region r = mesh_cube.newHex(n1, n2, n3, n4, n5, n6, n7, n8);
        
    }
    VTKWriter<IGMesh> writer_cube(mesh_cube);
    
    std::stringstream file_name_cube;
    file_name_cube<<"OLD_HLBFGS_DEBUG_CUBE_" << nb_file;
    writer_cube.write(file_name_cube.str(), DIM3 | R | N);
    
    
    nb_file++;

}
/*---------------------------------------------------------------------------*/

