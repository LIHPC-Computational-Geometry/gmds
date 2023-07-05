/*----------------------------------------------------------------------------*/
#include <fstream>
#include <bitset>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/PointSmoothing.h>
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/EdgeCollapse.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/EdgeInsertion.h>
#include <gmds/hybridMeshAdapt/DelaunayPointInsertion.h>
#include <gmds/hybridMeshAdapt/Octree.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  std::string fIn, fOut1, fOut2, fOut3, fOut4;
  if(argc != 2)
  {
      throw gmds::GMDSException("NO INPUT FILE : <mesh_file>");
  }
  fIn = std::string(argv[1]);
  std::string extansion(".vtk");
  std::size_t position = fIn.find(extansion);
  fOut1 = fIn.substr(0,position) +  "_HISTOGRAMME_TET_EUCLIDEAN.txt";
  fOut2 = fIn.substr(0,position) +  "_HISTOGRAMME_TET_METRIC.txt";
  fOut3 = fIn.substr(0,position) +  "_HISTOGRAMME_HEX.txt";
  fOut4 = fIn.substr(0,position) +  "_HISTOGRAMME_LENGTH.txt";
  if (fIn.find('.vtk') == std::string::npos) {
    throw gmds::GMDSException("<mesh_file> NOT A .vtk FILE");
  }
  std::cout << "INPUT FILE: " << fIn << std::endl;


  //==================================================================
  // MESH FILE READING
  //==================================================================
  Mesh m(MeshModel(DIM3 | R | F | E | N |
    R2N | F2N | E2N | R2F | F2R |
    F2E | E2F | R2E | N2R | N2F | N2E));

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N);
    vtkReader.read(fIn);

    Variable<math::Vector3d>* NODE_METRIC = nullptr;

    try{
      NODE_METRIC = m.getVariable<math::Vector3d, GMDS_NODE>("METRIC_NODES_HIST");
    }catch (gmds::GMDSException e)
    {
      //useful when no metric is specified or the metric is isotropic but not computed
      std::cout << "no metric attached" << std::endl;
      math::Vector3d mat({1.0/(1.36472*1.36472), 1.0/(1.36472*1.36472), 1.0/(1.36472*1.36472)});
      NODE_METRIC = m.newVariable<math::Vector3d, GMDS_NODE>("METRIC_NODES_HIST");
      NODE_METRIC->setValuesTo(mat);
      //throw gmds::GMDSException(e);
    }

    std::cout << " m.getNbTetrahedra() -> " << m.getNbTetrahedra() << std::endl;
    std::cout << " m.getNbHexahedra() -> " << m.getNbHexahedra() << std::endl;

  std::vector<double> hist_tet_euclidean{};
  std::vector<double> hist_tet_metric{};
  std::vector<double> hist_hex{};
  std::vector<double> hist_edge_length{};

  std::vector<std::pair<TInt, TInt>> edges{};
  std::set<std::pair<TInt, TInt>> seen{};
  auto make_pair = [&](const TInt indiceA, const TInt indiceB){
    return std::pair<TInt, TInt>{std::min(indiceA, indiceB), std::max(indiceA, indiceB)};
  };

  for(auto t:m.regions()){
    Region r = m.get<Region>(t);
    std::vector<TCellID> node_ids=r.getIDs<Node>();
    if(node_ids.size() == 8)
    {
      double qi = r.computeScaledJacobian();
      hist_hex.push_back(qi);
      std::pair<TInt, TInt> e01 = make_pair(node_ids[0], node_ids[1]);
      std::pair<TInt, TInt> e12 = make_pair(node_ids[1], node_ids[2]);
      std::pair<TInt, TInt> e23 = make_pair(node_ids[2], node_ids[3]);
      std::pair<TInt, TInt> e30 = make_pair(node_ids[3], node_ids[0]);

      std::pair<TInt, TInt> e45 = make_pair(node_ids[4], node_ids[5]);
      std::pair<TInt, TInt> e56 = make_pair(node_ids[5], node_ids[6]);
      std::pair<TInt, TInt> e67 = make_pair(node_ids[6], node_ids[7]);
      std::pair<TInt, TInt> e74 = make_pair(node_ids[7], node_ids[4]);

      std::pair<TInt, TInt> e04 = make_pair(node_ids[0], node_ids[4]);
      std::pair<TInt, TInt> e15 = make_pair(node_ids[1], node_ids[5]);
      std::pair<TInt, TInt> e26 = make_pair(node_ids[2], node_ids[6]);
      std::pair<TInt, TInt> e37 = make_pair(node_ids[3], node_ids[7]);

      if(seen.find(e01) == seen.end())
      edges.push_back(e01);
      if(seen.find(e12) == seen.end())
      edges.push_back(e12);
      if(seen.find(e23) == seen.end())
      edges.push_back(e23);
      if(seen.find(e30) == seen.end())
      edges.push_back(e30);

      if(seen.find(e45) == seen.end())
      edges.push_back(e45);
      if(seen.find(e56) == seen.end())
      edges.push_back(e56);
      if(seen.find(e67) == seen.end())
      edges.push_back(e67);
      if(seen.find(e74) == seen.end())
      edges.push_back(e74);

      if(seen.find(e04) == seen.end())
      edges.push_back(e04);
      if(seen.find(e15) == seen.end())
      edges.push_back(e15);
      if(seen.find(e26) == seen.end())
      edges.push_back(e26);
      if(seen.find(e37) == seen.end())
      edges.push_back(e37);

    }
    else if(node_ids.size() == 4)
    {
      Eigen::Matrix3d metric0;
      Eigen::Matrix3d metric1;
      Eigen::Matrix3d metric2;
      Eigen::Matrix3d metric3;
      metric0 << (*NODE_METRIC)[node_ids[0]].X(), 0.0 , 0.0 ,
      0.0 , (*NODE_METRIC)[node_ids[0]].Y() , 0.0 ,
      0.0 , 0.0 , (*NODE_METRIC)[node_ids[0]].Z();

      metric1 << (*NODE_METRIC)[node_ids[1]].X() , 0.0 , 0.0 ,
      0.0 , (*NODE_METRIC)[node_ids[1]].Y() , 0.0 ,
      0.0 , 0.0 , (*NODE_METRIC)[node_ids[1]].Z();

      metric2 << (*NODE_METRIC)[node_ids[2]].X() , 0.0 , 0.0 ,
      0.0 , (*NODE_METRIC)[node_ids[2]].Y() , 0.0 ,
      0.0 , 0.0 , (*NODE_METRIC)[node_ids[2]].Z();

      metric3 << (*NODE_METRIC)[node_ids[3]].X() , 0.0 , 0.0 ,
      0.0 , (*NODE_METRIC)[node_ids[3]].Y() , 0.0 ,
      0.0 , 0.0 , (*NODE_METRIC)[node_ids[3]].Z();

      double qi_m = r.computeQualityWithMetric(metric0, metric1,
                                                  metric2, metric3);
      double qi = r.computeQuality();
      hist_tet_euclidean.push_back(qi);
      hist_tet_metric.push_back(qi_m);

      std::pair<TInt, TInt> e01 = make_pair(node_ids[0], node_ids[1]);
      std::pair<TInt, TInt> e12 = make_pair(node_ids[1], node_ids[2]);
      std::pair<TInt, TInt> e20 = make_pair(node_ids[2], node_ids[0]);
      std::pair<TInt, TInt> e30 = make_pair(node_ids[3], node_ids[0]);
      std::pair<TInt, TInt> e31 = make_pair(node_ids[3], node_ids[1]);
      std::pair<TInt, TInt> e32 = make_pair(node_ids[3], node_ids[2]);


      if(seen.find(e01) == seen.end())
      edges.push_back(e01);
      if(seen.find(e12) == seen.end())
      edges.push_back(e12);
      if(seen.find(e20) == seen.end())
      edges.push_back(e20);
      if(seen.find(e30) == seen.end())
      edges.push_back(e30);
      if(seen.find(e31) == seen.end())
      edges.push_back(e31);
      if(seen.find(e32) == seen.end())
      edges.push_back(e32);
    }
  }

  for(auto const & edge : edges){
    Node n0 = m.get<Node>(edge.first);
    Node n1 = m.get<Node>(edge.second);
    math::Point p0 = n0.point();
    math::Point p1 = n1.point();
    math::Vector3d v01_ =  p0 - p1;
    Eigen::Vector3d v01 = Eigen::Vector3d(v01_.X(), v01_.Y(), v01_.Z());

    Eigen::Matrix3d metric0;
    Eigen::Matrix3d metric1;
    metric0 << (*NODE_METRIC)[edge.first].X(), 0.0 , 0.0 ,
    0.0 , (*NODE_METRIC)[edge.first].Y() , 0.0 ,
    0.0 , 0.0 , (*NODE_METRIC)[edge.first].Z();

    metric1 << (*NODE_METRIC)[edge.second].X() , 0.0 , 0.0 ,
    0.0 , (*NODE_METRIC)[edge.second].Y() , 0.0 ,
    0.0 , 0.0 , (*NODE_METRIC)[edge.second].Z();

    double length_metric = 0.5*sqrt(v01.transpose() * (metric0*v01)) + 0.5*sqrt(v01.transpose() * (metric1*v01));
    hist_edge_length.push_back(length_metric);
  }
  std::cout << " hist_tet_euclidean -> " << hist_tet_euclidean.size() << std::endl;
  std::cout << " hist_tet_metric -> " << hist_tet_metric.size() << std::endl;
  std::cout << " hist_hex -> " << hist_hex.size() << std::endl;
  std::cout << " hist_edge_length -> " << hist_edge_length.size() << std::endl;
  // Création des objets ofstream
  std::ofstream fichier_tet_euclidean;
  std::ofstream fichier_tet_metric;
  std::ofstream fichier_hex;
  std::ofstream fichier_length_edges;

  fichier_tet_euclidean.open(fOut1, std::ios::out);
  fichier_tet_metric.open(fOut2, std::ios::out);
  fichier_hex.open(fOut3, std::ios::out);
  fichier_length_edges.open(fOut4, std::ios::out);

  // Si erreur d'ouverture
  if(fichier_tet_euclidean.bad() || fichier_tet_metric.bad() || fichier_hex.bad())
      return 0; // on quitte

  for(auto const v : hist_tet_euclidean)
    fichier_tet_euclidean << v << std::endl;

  for(auto const v : hist_tet_metric)
      fichier_tet_metric << v << std::endl;

  for(auto const v : hist_hex)
    fichier_hex << v << std::endl;

  for(auto const v : hist_edge_length)
    fichier_length_edges << v << std::endl;

  // Fermeture des fichier
  fichier_tet_euclidean.close();
  fichier_tet_metric.close();
  fichier_hex.close();
  fichier_length_edges.close();

  return 0;
}

/*----------------------------------------------------------------------------*/
