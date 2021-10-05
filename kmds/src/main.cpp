
#include <Kokkos_Core.hpp>
#include <cstdio>
#include <typeinfo>

#include <KM/DS/Connectivity.h>
#include <KM/DS/Mesh.h>
#include <KM/Utils/GrowingView.h>

int nb_nodes_i;
int nb_nodes_j;
int nb_faces_i;
int nb_faces_j;
int nb_nodes;
int nb_faces;


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

struct SelectFace
{
        kmds::Mesh* m;
        Filter* filter;

        kmds::GrowingView<kmds::TCellID>* result;
        // Constructor takes View by "value"; this does a shallow copy.
        SelectFace(kmds::Mesh* m_, Filter* f_, kmds::GrowingView<kmds::TCellID>* r_)
         : m(m_)
         , filter(f_)
         , result(r_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
                kmds::Face f = m->getFace(i);
                if (filter->select(f.id)) {
                        //    printf("%i) %i is selected!!!\n", i, f.id);
                        result->push_back(f.id);
                }
        }
};

struct Q2T
{
        kmds::Mesh* m;

        kmds::GrowingView<kmds::TCellID>* faces;
        // Constructor takes View by "value"; this does a shallow copy.
        Q2T(kmds::Mesh* m_, kmds::GrowingView<kmds::TCellID>* f_)
         : m(m_)
         , faces(f_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
                int fid = faces->get(i);
                // printf("> Face %i\n", fid);
                kmds::Face f = m->getFace(fid);
                // Kokkos::View<kmds::Node*> n("nf", f.getNbNodes());
                // f.nodes(n);
//                Kokkos::View<kmds::Node*> n = f.nodes();
//
//                // We reuse face node spaces to put the first new triangle
//
//                Kokkos::View<kmds::TCellID*> t1("T1", 3);
//                t1(0) = n(0).id;
//                t1(1) = n(1).id;
//                t1(2) = n(2).id;
//
//                f.setNodes(t1);
//                // and we add the second one
//                m->newTriangle(n(2).id, n(3).id, n(0).id);
                kmds::Node n[5];
                int nsize = 0;
                f.nodes(n, &nsize);

                // We reuse face node spaces to put the first new triangle
                kmds::TCellID t1[5];
                t1[0] = n[0].id;
                t1[1] = n[1].id;
                t1[2] = n[2].id;

                f.setNodes(t1,3);

                // and we add the second one
                m->newTriangle(n[2].id, n[3].id, n[0].id);
        }
};

//struct ColorFace
//{
//        kmds::Mesh* m;
//
//        kmds::Variable<int>* v;
//        // Constructor takes View by "value"; this does a shallow copy.
//        ColorFace(kmds::Mesh* m_, kmds::Variable<int>* v_)
//         : m(m_)
//         , v(v_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Face f = m->getFace(i);
//                if (f.id % 2 == 0) {
//                        (*v)[f] = 2;
//                } else {
//                        (*v)[f] = 0;
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
//                kmds::Face f = m->getFace(i);
//                printf("Color of %i: %i\n", i, (*v)[f]);
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
                double x = (a % 4);
                double y = (a / 4);
                kmds::Node n = m->getNode(a);
                n.setLocation(x, y, 0);
                //   printf("For %i, create node %i at (%f,%f)\n", i, a, x, y);
        }
};
struct CreateQuad
{
        kmds::Mesh* m;
        // Constructor takes View by "value"; this does a shallow copy.
        CreateQuad(kmds::Mesh* m_)
         : m(m_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
                int a = m->addQuad();
                int x = (a % 3) + 4 * (a / 3);
                int x2 = x + 1;
                int x3 = x + 4;
                int x4 = x + 5;
                kmds::Face q = m->getFace(a);
//                Kokkos::View<kmds::TCellID*> v("qn", 4);
//                v(0) = x;
//                v(1) = x2;
//                v(2) = x4;
//                v(3) = x3;
//                q.setUnsafeNodes(v);

                kmds::TCellID v[5];
                v[0] = x;
                v[1] = x2;
                v[2] = x4;
                v[3] = x3;
                q.setNodes(v,4);

                //    printf("For %i, create quad %i with (%i, %i, %i, %i)\n", i, a, x, x2, x3, x4);
        }
};

//struct GetQuad
//{
//        kmds::Mesh* m;
//        // Constructor takes View by "value"; this does a shallow copy.
//        GetQuad(kmds::Mesh* m_)
//         : m(m_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Face q = m->getFace(i);
//                // std::vector<kmds::TCellID> v;
//                // q.idNodes(v);
//                // printf("For %i, get quad %i with (%i, %i, %i, %i)\n", i, i, v[0], v[1], v[2], v[3]);
//
//                Kokkos::View<kmds::TCellID*> kv;
//                q.nodeIds(kv);
//                printf("For %i, get quad %i with (%i, %i, %i, %i)\n", i, i, kv(0), kv(1), kv(2), kv(3));
//        }
//};
//struct GetQuadByObject
//{
//        kmds::Mesh* m;
//        // Constructor takes View by "value"; this does a shallow copy.
//        GetQuadByObject(kmds::Mesh* m_)
//         : m(m_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Face q = m->getFace(i);
//                // std::vector<kmds::TCellID> v;
//                // q.idNodes(v);
//                // printf("For %i, get quad %i with (%i, %i, %i, %i)\n", i, i, v[0], v[1], v[2], v[3]);
//
//                Kokkos::View<kmds::Node*> kv("knodes", q.getNbNodes());
//                q.nodes(kv);
//                printf("For %i, get quad by object %i with (%i, %i, %i, %i)\n", i, i, kv(0).id, kv(1).id, kv(2).id,
//                       kv(3).id);
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
//                printf("%i) Point (%f,%f)\n", i, x, y);
//        }
//};
//
//struct DisplayN2F
//{
//        kmds::Mesh* m;
//        // Constructor takes View by "value"; this does a shallow copy.
//        DisplayN2F(kmds::Mesh* m_)
//         : m(m_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const
//        {
//                kmds::Node n = m->getNode(i);
//                Kokkos::View<kmds::TCellID*> fids;
//                n.faceIds(fids);
//                int k = fids.size() + 1;
//                printf("Node %i has %i adj faces\n", n.id, k);
//                for (auto i = 0; i <= fids.size(); i++) {
//                        int j = fids(i);
//                        printf("N %i --> F %i\n", n.id, j);
//                }
//        }
//};

int
main(int argc, char* argv[])
{
        int num_threads;
        int use_numa = -1;
        int use_core = -1;

        if (argc == 4) {
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

        } else {
                std::cerr<<"need number of threads argument."<<std::endl;
                exit(-1);
        }

        nb_faces_i = nb_nodes_i - 1;
        nb_faces_j = nb_nodes_j - 1;
        nb_nodes = nb_nodes_i * nb_nodes_j;
        nb_faces = nb_faces_i * nb_faces_j;

//        int nb_nodes = 1000000;
//        int nb_faces = 998001;
        std::cout << "Kokkos::init" << std::endl;
	Kokkos::InitArguments kargs;
	kargs.num_threads = num_threads;
	//kargs.use_numa = use_numa;
	//kargs.use_core = use_core;
	Kokkos::initialize(kargs);
	
      
        Kokkos::Timer timer;
        kmds::Mesh m;
        m.updateFaceCapacity(200000000);
        m.updateNodeCapacity(200000000);
        // //===========================================================
        // Create a regular grid of nodes (in 2D)
        Kokkos::parallel_for(nb_nodes, AddNode(&m));

        std::cout << "Nodes (nb/C): (" << m.getNbNodes() << ", " << m.getNodeCapacity() << ")" << std::endl;

        //===========================================================
        // Display all mesh nodes
        // Kokkos::parallel_for(100, DisplayNode(&m));

        //===========================================================
        // Create quads in the grid
        Kokkos::parallel_for(nb_faces, CreateQuad(&m));

        //===========================================================
        // Access to all the quads (id + id vertex list) as ids
        // Kokkos::parallel_for(nb_faces, GetQuad(&m));
        //===========================================================
        // Access to all the quads (id + id vertex list) as Face
        // handler
        //  Kokkos::parallel_for(9801, GetQuadByObject(&m));

        timer.reset();

        //===========================================================
        // We store in a view faces that verify a certain criteria
        printf("==== NB FACES = %i ======\n", (int) m.getNbFaces());

        kmds::GrowingView<kmds::TCellID> selection("SELECTION", m.getNbFaces());
        Filter f(true);
        Kokkos::parallel_for(m.getNbFaces(), SelectFace(&m, &f, &selection));
        printf("Nb selected faces = %i\n", selection.getNbElems());

        //===========================================================
        // Each stored face (quad) is now split into 2 triangles
        Kokkos::parallel_for(selection.getNbElems(), Q2T(&m, &selection));
        int nbf = m.getNbFaces();
//        int nbq = m.getNbQuads();
//        int nbt = m.getNbTriangles();
//        printf("(Faces, Quads, Tris) = (%i,%i,%i)\n", nbf, nbq, nbt);
        printf("(Faces) = (%i)\n", nbf);

        //===========================================================
        // Now we create a variable that lives on faces and we use it
//        kmds::Variable<int>* v = m.createVariable<int>(kmds::KMDS_FACE, "color");
//        Kokkos::parallel_for(m.getNbFaces(), ColorFace(&m, v));
        // Kokkos::parallel_for(m.getNbFaces(), DisplayColor(&m));
    printf("nbThreads_Time: %i %f\n", num_threads, timer.seconds());

        // kmds::Connectivity* N2F = m.createConnectivity(kmds::N2F);
        // kmds::ConnectivityHelper ch(&m);
        // ch.buildN2F();
        // Kokkos::parallel_for(m.getNbNodes(), DisplayN2F(&m));

       Kokkos::finalize();
}
