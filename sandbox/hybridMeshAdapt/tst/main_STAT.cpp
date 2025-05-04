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
  std::string fIn, fOut;
  if(argc != 2)
  {
      throw gmds::GMDSException("NO INPUT FILE : <mesh_file> ");
  }
  fIn = std::string(argv[1]);
  fOut ="carateristics_cells.txt";
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
    vtkReader.read(fIn);

  //compute the hex number
  const unsigned int tetNbr = m.getNbTetrahedra();
  const unsigned int tetHexa = m.getNbHexahedra();

  double vol_tet = 0.0;
  double vol_hex = 0.0;
  int hexCpt = 0;
  int tetCpt = 0;
  for(auto t:m.regions()){
    Region r = m.get<Region>(t);
    std::vector<TCellID> node_ids=r.getIDs<Node>();
    if(node_ids.size() == 8)
    {
      vol_hex += std::abs(r.volume());
    }
    else if(node_ids.size() == 4)
    {
      vol_tet += std::abs(r.volume());
    }
  }

  const double voltot = vol_hex + vol_tet;
  const double percVolTet = vol_tet / voltot;
  const double percVolHex = vol_hex / voltot;

  const unsigned int tot = tetNbr + tetHexa;
  const double percentageTet = static_cast<double>(tetNbr) / static_cast<double>(tot);
  const double percentageHex = static_cast<double>(tetHexa) / static_cast<double>(tot);

  std::cout << "tetNbr -> " << tetNbr << std::endl;
  std::cout << "tetHexa -> " << tetHexa << std::endl;

  std::cout << "percentageTet -> " << percentageTet << std::endl;
  std::cout << "percentageHex -> " << percentageHex << std::endl;

  std::cout << "percVolTet -> " << percVolTet << std::endl;
  std::cout << "percVolHex -> " << percVolHex << std::endl;

  // Création d'un objet ofstream -> pour ecrire
  std::ofstream fichier;
  fichier.open(fOut, std::ios::out  | std::ios::app);

  // Si erreur d'ouverture
  if(fichier.bad())
      return 0; // on quitte

  fichier << "fIn -> " << fIn << std::endl;
  fichier << "tetNbr -> " << tetNbr << std::endl;
  fichier << "tetHexa -> " << tetHexa << std::endl;
  fichier << "percentageTet -> " << percentageTet << std::endl;
  fichier << "percentageHex -> " << percentageHex << std::endl;
  fichier << "percVolTet -> " << percVolTet << std::endl;
  fichier << "percVolHex -> " << percVolHex << std::endl;
  fichier << std::endl;
  fichier << std::endl;
  fichier << std::endl;
  fichier << std::endl;
  // Fermeture du fichier
  fichier.close();

  return 0;
}

/*----------------------------------------------------------------------------*/
