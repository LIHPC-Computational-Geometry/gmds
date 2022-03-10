
#include <Kokkos_Core.hpp>
#include <cstdio>
#include <typeinfo>

#include <KM/DS/Connectivity.h>
#include <KM/DS/Mesh.h>
#include <KM/Utils/GrowingView.h>

//const int nb_nodes_i = 100;
//const int nb_nodes_j = 100;
//const int nb_nodes_k = 100;
//const int nb_regions_i = nb_nodes_i - 1;
//const int nb_regions_j = nb_nodes_j - 1;
//const int nb_regions_k = nb_nodes_k - 1;
//const int nb_nodes = nb_nodes_i * nb_nodes_j * nb_nodes_k;
//const int nb_regions = nb_regions_i * nb_regions_j * nb_regions_k;
int nb_nodes_i;
int nb_nodes_j;
int nb_nodes_k;
int nb_regions_i;
int nb_regions_j;
int nb_regions_k;
int nb_nodes;
int nb_regions;

class Filter
{
 public:
        Filter(const bool AEven)
         : even(AEven)
        {
                ;
        }

        bool
        select(kmds::TCellID AID)
        {
                if (even && (AID % 2 == 0))
                        return true;
                else if (!even && (AID % 2 != 0))
                        return true;

                return false;
        }

 private:
        bool even;
};

struct SelectRegion
{
        kmds::Mesh* m;
        Filter* filter;

        kmds::GrowingView<kmds::TCellID>* result;
        // Constructor takes View by "value"; this does a shallow copy.
        SelectRegion(kmds::Mesh* m_, Filter* f_, kmds::GrowingView<kmds::TCellID>* r_)
         : m(m_)
         , filter(f_)
         , result(r_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
                kmds::Region r = m->getRegion(i);
                if (filter->select(r.id)) {
                        //    printf("%i) %i is selected!!!\n", i, r.id);
                        result->push_back(r.id);
                }
        }
};

struct H2P
{
        kmds::Mesh* m;

        kmds::GrowingView<kmds::TCellID>* regions;
        // Constructor takes View by "value"; this does a shallow copy.
        H2P(kmds::Mesh* m_, kmds::GrowingView<kmds::TCellID>* f_)
         : m(m_)
         , regions(f_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
                int rid = regions->get(i);
                //printf("> Region %i\n", rid);
                kmds::Region r = m->getRegion(rid);
                // Kokkos::View<kmds::Node*> n("nr", r.getNbNodes());
                // r.nodes(n);
                //Kokkos::View<kmds::Node*> n = r.nodes();
                kmds::Node n[12];
                int nsize = 0;
                r.nodes(n, &nsize);

                // we create a node
                kmds::TCoord x = 0.;
                kmds::TCoord y = 0.;
                kmds::TCoord z = 0.;
                
                for(int i=0; i<8; i++) {
                        kmds::TCoord xtmp;
                        kmds::TCoord ytmp;
                        kmds::TCoord ztmp;
                        n[i].getLocation(xtmp, ytmp, ztmp);
                        x += xtmp;
                        y += ytmp;
                        z += ztmp;
                }
                
                kmds::TCellID newn = m->newNode(x/8., y/8., z/8.);

                // We reuse region node spaces to put the first new pyramid

//                Kokkos::View<kmds::TCellID*> p1("P1", 5); // bottom
//                p1(0) = n(0).id;
//                p1(1) = n(1).id;
//                p1(2) = n(2).id;
//                p1(3) = n(3).id;
//                p1(4) = newn;
//                r.setNodes(p1);
//                // and we add the five other ones
//                m->newPyramid(n(4).id, n(7).id, n(6).id, n(5).id, newn); // top
//                m->newPyramid(n(0).id, n(3).id, n(7).id, n(4).id, newn); // left
//                m->newPyramid(n(2).id, n(1).id, n(5).id, n(6).id, newn); // right
//                m->newPyramid(n(1).id, n(0).id, n(4).id, n(5).id, newn); // front
//                m->newPyramid(n(3).id, n(2).id, n(6).id, n(7).id, newn); // back

                kmds::TCellID p1[5];
                p1[0] = n[0].id;
                p1[1] = n[1].id;
                p1[2] = n[2].id;
                p1[3] = n[3].id;
                p1[4] = newn;
                r.setNodes(p1,5);
                // and we add the five other ones
                m->newPyramid(n[4].id, n[7].id, n[6].id, n[5].id, newn); // top
                m->newPyramid(n[0].id, n[3].id, n[7].id, n[4].id, newn); // left
                m->newPyramid(n[2].id, n[1].id, n[5].id, n[6].id, newn); // right
                m->newPyramid(n[1].id, n[0].id, n[4].id, n[5].id, newn); // front
                m->newPyramid(n[3].id, n[2].id, n[6].id, n[7].id, newn); // back
        }
};

//struct ColorRegion
//{
//        kmds::Mesh* m;
//
//        kmds::Variable<int>* v;
//        // Constructor takes View by "value"; this does a shallow copy.
//        ColorRegion(kmds::Mesh* m_, kmds::Variable<int>* v_)
//         : m(m_)
//         , v(v_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Region r = m->getRegion(i);
//                if (r.id % 2 == 0) {
//                        (*v)[r] = 2;
//                } else {
//                        (*v)[r] = 0;
//                }
//        }
//};
//struct DisplayColor
//{
//        kmds::Mesh* m;
//        kmds::Variable<int>* v;
//        // Constructor takes View by "value"; this does a shallow copy.
//        DisplayColor(kmds::Mesh* m_)
//         : m(m_)
//        {
//                v = m->getVariable<int>(kmds::KMDS_NODE, "color");
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Region r = m->getRegion(i);
//                printf("Color of %i: %i\n", i, (*v)[r]);
//        }
//};

struct AddNode
{
        kmds::Mesh* m;
        // Constructor takes View by "value"; this does a shallow copy.
        AddNode(kmds::Mesh* m_)
         : m(m_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
                int a = m->addNode();
                kmds::Node n = m->getNode(a);
                int x = floor((double) a / ((double) (nb_nodes_j * nb_nodes_k)));
                a -= x * (nb_nodes_j * nb_nodes_k);
                int y = floor((double) a / (double) nb_nodes_k);
                a -= y * nb_nodes_k;
                int z = a;
                n.setLocation(x, y, z);
                //printf("For %i, create node %i at (%f,%f,%f)\n", i, a, (double)x, (double)y, (double)z);
        }
};
struct CreateHex
{
        kmds::Mesh* m;
        // Constructor takes View by "value"; this does a shallow copy.
        CreateHex(kmds::Mesh* m_)
         : m(m_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
                int index = m->addHexahedron();
                int a = index;
                int x = floor((double) a / ((double) (nb_regions_j * nb_regions_k)));
                a -= x * (nb_regions_j * nb_regions_k);
                int y = floor((double) a / (double) nb_regions_k);
                a -= y * nb_regions_k;
                int z = a;
 
                int x1 = x     * (nb_nodes_j * nb_nodes_k) + y     * nb_nodes_k + z;
                int x2 = (x+1) * (nb_nodes_j * nb_nodes_k) + y     * nb_nodes_k + z;
                int x3 = (x+1) * (nb_nodes_j * nb_nodes_k) + (y+1) * nb_nodes_k + z;
                int x4 = x     * (nb_nodes_j * nb_nodes_k) + (y+1) * nb_nodes_k + z;
                int x5 = x1 + 1;
                int x6 = x2 + 1;
                int x7 = x3 + 1;
                int x8 = x4 + 1;
                kmds::Region r = m->getRegion(a);
//                Kokkos::View<kmds::TCellID*> v("qn", 8);
//                v(0) = x1;
//                v(1) = x2;
//                v(2) = x4;
//                v(3) = x3;
//                v(4) = x5;
//                v(5) = x6;
//                v(6) = x7;
//                v(7) = x8;
//                r.setUnsafeNodes(v);
                kmds::TCellID v[12];
                v[0] = x1;
                v[1] = x2;
                v[2] = x4;
                v[3] = x3;
                v[4] = x5;
                v[5] = x6;
                v[6] = x7;
                v[7] = x8;
                r.setNodes(v,8);
                //printf("For %i, create hex %i with (%i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i)\n", i, index, x,y,z,x1, x2, x3, x4, x5, x6, x7, x8);
        }
};

//struct GetHex
//{
//        kmds::Mesh* m;
//        // Constructor takes View by "value"; this does a shallow copy.
//        GetHex(kmds::Mesh* m_)
//         : m(m_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Region r = m->getRegion(i);
//                // std::vector<kmds::TCellID> v;
//                // q.idNodes(v);
//                // printf("For %i, get quad %i with (%i, %i, %i, %i)\n", i, i, v[0], v[1], v[2], v[3]);
//
//                Kokkos::View<kmds::TCellID*> kv;
//                r.nodeIds(kv);
//                printf("For %i, get hex %i with (%i, %i, %i, %i, %i, %i, %i, %i)\n", i, i, kv(0), kv(1), kv(2), kv(3), kv(4), kv(5), kv(6), kv(7));
//        }
//};
//struct GetHexByObject
//{
//        kmds::Mesh* m;
//        // Constructor takes View by "value"; this does a shallow copy.
//        GetHexByObject(kmds::Mesh* m_)
//         : m(m_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Region r = m->getRegion(i);
//                // std::vector<kmds::TCellID> v;
//                // q.idNodes(v);
//                // printf("For %i, get quad %i with (%i, %i, %i, %i)\n", i, i, v[0], v[1], v[2], v[3]);
//
//                Kokkos::View<kmds::Node*> kv("knodes", r.getNbNodes());
//                r.nodes(kv);
//                printf("For %i, get hex by object %i with (%i, %i, %i, %i, %i, %i, %i, %i)\n", i, i, kv(0).id, kv(1).id, kv(2).id,
//                       kv(3).id, kv(4).id, kv(5).id, kv(6).id, kv(7).id);
//        }
//};
//
//struct DisplayNode
//{
//        kmds::Mesh* m;
//        // Constructor takes View by "value"; this does a shallow copy.
//        DisplayNode(kmds::Mesh* m_)
//         : m(m_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Node n = m->getNode(i);
//                double x, y, z;
//                n.getLocation(x, y, z);
//                printf("%i) Point (%f,%f,%f)\n", i, x, y, z);
//        }
//};
//
//struct DisplayN2R
//{
//        kmds::Mesh* m;
//        // Constructor takes View by "value"; this does a shallow copy.
//        DisplayN2R(kmds::Mesh* m_)
//         : m(m_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Node n = m->getNode(i);
//                Kokkos::View<kmds::TCellID*> rids;
//                n.regionIds(rids);
//                int k = rids.size() + 1;
//                printf("Node %i has %i adj regions\n", n.id, k);
//                for (auto i = 0; i <= rids.size(); i++) {
//                        int j = rids(i);
//                        printf("N %i --> R %i\n", n.id, j);
//                }
//        }
//};

int
main(int argc, char* argv[])
{
        int num_threads;
        int use_numa = -1;
        int use_core = -1;

        if (argc == 5) {
                std::istringstream iss( argv[1] );
                int val;

                if (iss >> val)
                {
                        num_threads = val;
                } else {
                        std::cerr<<"could not convert number of threads argument."<<std::endl;
                        exit(-1);
                }

            std::istringstream iss_2( argv[2] );
            iss_2 >> nb_nodes_i;
            std::istringstream iss_3( argv[3] );
            iss_3 >> nb_nodes_j;
            std::istringstream iss_4( argv[4] );
            iss_4 >> nb_nodes_k;

        } else {
                std::cerr<<"need nb_threads nb_nodes_i nb_nodes_j arguments."<<std::endl;
                exit(-1);
        }

    nb_regions_i = nb_nodes_i - 1;
    nb_regions_j = nb_nodes_j - 1;
    nb_regions_k = nb_nodes_k - 1;
    nb_nodes = nb_nodes_i * nb_nodes_j * nb_nodes_k;
    nb_regions = nb_regions_i * nb_regions_j * nb_regions_k;


      	std::cout << "Kokkos::init" << std::endl;
	Kokkos::InitArguments kargs;
	kargs.num_threads = num_threads;
	//kargs.use_numa = use_numa;
	//kargs.use_core = use_core;
	Kokkos::initialize(kargs);

	
        Kokkos::Timer timer;
        kmds::Mesh m;
        m.updateRegionCapacity(100000000);
        m.updateNodeCapacity(20000000);
        // //===========================================================
        // Create a regular grid of nodes (in 3D)
        Kokkos::parallel_for(nb_nodes, AddNode(&m));

        std::cout << "Nodes (nb/C): (" << m.getNbNodes() << ", " << m.getNodeCapacity() << ")" << std::endl;

        //===========================================================
        // Display all mesh nodes
        // Kokkos::parallel_for(100, DisplayNode(&m));

        //===========================================================
        // Create quads in the grid
        Kokkos::parallel_for(nb_regions, CreateHex(&m));

        //===========================================================
        // Access to all the quads (id + id vertex list) as ids
        // Kokkos::parallel_for(nb_faces, GetQuad(&m));
        //===========================================================
        // Access to all the quads (id + id vertex list) as Face
        // handler
        //  Kokkos::parallel_for(9801, GetQuadByObject(&m));

        timer.reset();

        //===========================================================
        // We store in a view regions that verify a certain criteria
        printf("==== NB REGIONS = %i ======\n", (int) m.getNbRegions());
        kmds::GrowingView<kmds::TCellID> selection("SELECTION", m.getNbRegions());
        Filter f(true);
        Kokkos::parallel_for(m.getNbRegions(), SelectRegion(&m, &f, &selection));
        printf("Nb selected regions = %i\n", selection.getNbElems());

        //===========================================================
        // Each stored region (hex) is now split into 6 pyramids
        Kokkos::parallel_for(selection.getNbElems(), H2P(&m, &selection));
        int nbr = m.getNbRegions();
//        int nbh = m.getNbHexahedra();
//        int nbp = m.getNbPyramids();
//        printf("(Regions, Hexahedra, Pyramids) = (%i,%i,%i)\n", nbr, nbh, nbp);
        printf("(Regions) = (%i)\n", nbr);
        //===========================================================
        // Now we create a variable that lives on faces and we use it
//        kmds::Variable<int>* v = m.createVariable<int>(kmds::KMDS_REGION, "color");
//        Kokkos::parallel_for(m.getNbRegions(), ColorRegion(&m, v));
        // Kokkos::parallel_for(m.getNbFaces(), DisplayColor(&m));
    printf("nbThreads_Time: %i %f\n", num_threads, timer.seconds());

        // kmds::Connectivity* N2F = m.createConnectivity(kmds::N2F);
        // kmds::ConnectivityHelper ch(&m);
        // ch.buildN2F();
        // Kokkos::parallel_for(m.getNbNodes(), DisplayN2F(&m));

    Kokkos::finalize();
}
