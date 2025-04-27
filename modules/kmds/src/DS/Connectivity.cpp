/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Connectivity.cpp
 *  \author  F. LEDOUX
 *  \date    03/13/2017
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// STL headers
#include <set>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <KM/DS/CellHandle.h>
#include <KM/DS/Connectivity.h>
#include <KM/DS/Mesh.h>
#include <Kokkos_UnorderedMap.hpp>
#include <stdlib.h>
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*------------------------------------------------------------------------*/
struct BuildN2EFunctor
{
    kmds::GrowingView<kmds::TCellID>* cellIDsAccessor;
    ConnectivityHelper* chelper;
    Mesh* m;
    Kokkos::View<TCellID * [10]> n2e;
    Kokkos::View<int*> n2nb;

    /*------------------------------------------------------------------------*/
    /** \brief Constructor
     *
     * \param[in] ACH the connectivty helper that calls this functor
     * \param[in] AM  the mesh it works on
     */
    BuildN2EFunctor(kmds::GrowingView<kmds::TCellID>* cellIDsAccessor_, ConnectivityHelper* ACH, Mesh* AM, Kokkos::View<TCellID * [10]> AN2E, Kokkos::View<int*> AN2NB)
            : cellIDsAccessor(cellIDsAccessor_)
            , chelper(ACH)
            , m(AM)
            , n2e(AN2E)
            , n2nb(AN2NB)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(TCellID i) const
    {
            Edge e = m->getEdge(cellIDsAccessor->get(i));
            Kokkos::View<TCellID*> n;
            e.nodeIds(n);
            for (auto i_n = 0; i_n < n.size(); i_n++) {
                    TCellID n_id = n(i_n);
                    int top = Kokkos::atomic_fetch_add(&n2nb(n_id), 1);
                    n2e(n_id, top) = e.id;
            }
    }
};
/*------------------------------------------------------------------------*/
struct BuildN2FFunctor
{
        kmds::GrowingView<kmds::TCellID>* cellIDsAccessor;
        ConnectivityHelper* chelper;
        Mesh* m;
        Kokkos::View<TCellID * [30]> n2f;
        Kokkos::View<int*> n2nb;

        /*------------------------------------------------------------------------*/
        /** \brief Constructor
         *
         * \param[in] ACH the connectivty helper that calls this functor
         * \param[in] AM  the mesh it works on
         */
        BuildN2FFunctor(kmds::GrowingView<kmds::TCellID>* cellIDsAccessor_, ConnectivityHelper* ACH, Mesh* AM, Kokkos::View<TCellID * [30]> AN2F, Kokkos::View<int*> AN2NB)
         : cellIDsAccessor(cellIDsAccessor_)
         , chelper(ACH)
         , m(AM)
         , n2f(AN2F)
         , n2nb(AN2NB)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(TCellID i) const
        {
                Face f = m->getFace(cellIDsAccessor->get(i));
                Kokkos::View<TCellID*> n;
                f.nodeIds(n);
                for (auto i_n = 0; i_n < n.size(); i_n++) {
                        TCellID n_id = n(i_n);
                        int top = Kokkos::atomic_fetch_add(&n2nb(n_id), 1);
                        n2f(n_id, top) = f.id;
                }
        }
};
/*------------------------------------------------------------------------*/
struct BuildN2RFunctor
{
    kmds::GrowingView<kmds::TCellID>* cellIDsAccessor;
    ConnectivityHelper* chelper;
    Mesh* m;
    Kokkos::View<TCellID * [20]> n2c;
    Kokkos::View<int*> n2nb;

    /*------------------------------------------------------------------------*/
    /** \brief Constructor
     *
     * \param[in] ACH the connectivty helper that calls this functor
     * \param[in] AM  the mesh it works on
     */
    BuildN2RFunctor(kmds::GrowingView<kmds::TCellID>* cellIDsAccessor_, ConnectivityHelper* ACH, Mesh* AM, Kokkos::View<TCellID * [20]> AN2C, Kokkos::View<int*> AN2NB)
            : cellIDsAccessor(cellIDsAccessor_)
            , chelper(ACH)
            , m(AM)
            , n2c(AN2C)
            , n2nb(AN2NB)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(TCellID i) const
    {
            Region r = m->getRegion(cellIDsAccessor->get(i));
            Kokkos::View<TCellID*> n;
            r.nodeIds(n);
            for (auto i_n = 0; i_n < n.size(); i_n++) {
                    TCellID n_id = n(i_n);
                    int top = Kokkos::atomic_fetch_add(&n2nb(n_id), 1);
                    n2c(n_id, top) = r.id;
            }
    }


};
/*------------------------------------------------------------------------*/
struct BuildN2RFunctor_variant_0
{
    ConnectivityHelper* chelper;
    Mesh* m;
    Kokkos::View<TCellID * [20]> n2c;
    Kokkos::View<int*> n2nb;

    /*------------------------------------------------------------------------*/
    /** \brief Constructor
     *
     * \param[in] ACH the connectivty helper that calls this functor
     * \param[in] AM  the mesh it works on
     */
    BuildN2RFunctor_variant_0(ConnectivityHelper* ACH, Mesh* AM, Kokkos::View<TCellID * [20]> AN2C, Kokkos::View<int*> AN2NB)
            : chelper(ACH)
            , m(AM)
            , n2c(AN2C)
            , n2nb(AN2NB)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(TCellID i) const
    {
            if(m->hasRegion(i)) {
                    Region r = m->getRegion(i);
                    Kokkos::View<TCellID *> n;
                    r.nodeIds(n);
                    for (auto i_n = 0; i_n < n.size(); i_n++) {
                            TCellID n_id = n(i_n);
                            int top = Kokkos::atomic_fetch_add(&n2nb(n_id), 1);
                            n2c(n_id, top) = r.id;
                    }
            }
    }


};
/*------------------------------------------------------------------------*/
struct BuildF2F_byNFunctor
{
    kmds::GrowingView<kmds::TCellID>* cellIDsAccessor;
    ConnectivityHelper* chelper;
    Mesh* m;
    Kokkos::View<TCellID * [20]> c2c;
    Kokkos::View<int*> c2nb;

    const kmds::Connectivity* c_N2F;

    /*------------------------------------------------------------------------*/
    /** \brief Constructor
     *
     * \param[in] ACH the connectivty helper that calls this functor
     * \param[in] AM  the mesh it works on
     */
    BuildF2F_byNFunctor(kmds::GrowingView<kmds::TCellID>* cellIDsAccessor_, ConnectivityHelper* ACH, Mesh* AM, const kmds::Connectivity* Ac_N2F, Kokkos::View<TCellID * [20]> AC2C, Kokkos::View<int*> AC2NB)
            : cellIDsAccessor(cellIDsAccessor_)
            , chelper(ACH)
            , m(AM)
            , c_N2F(Ac_N2F)
            , c2c(AC2C)
            , c2nb(AC2NB)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(TCellID i) const
    {
            TCellID r_id = cellIDsAccessor->get(i);
            Face f = m->getFace(r_id);
            Kokkos::View<TCellID*> n;
            f.nodeIds(n);

            std::set<TCellID> facesSet;

            for (auto i_n = 0; i_n < n.size(); i_n++) {
                    TCellID n_id = n(i_n);

                    Kokkos::View<TCellID*> faces_tmp;
                    c_N2F->get(n_id, faces_tmp);

                    for(int i_r=0; i_r<faces_tmp.size(); i_r++) {
                            facesSet.insert(faces_tmp(i_r));
                    }
            }
            facesSet.erase(r_id);

            c2nb(r_id)  = facesSet.size();

            int index = 0;
            for(auto r_tmp: facesSet) {
                    c2c(r_id, index) = r_tmp;
                    index++;
            }
    }
};
/*------------------------------------------------------------------------*/
struct BuildR2R_byNFunctor
{
    kmds::GrowingView<kmds::TCellID>* cellIDsAccessor;
    ConnectivityHelper* chelper;
    Mesh* m;
    Kokkos::View<TCellID * [40]> c2c;
    Kokkos::View<int*> c2nb;

    const kmds::Connectivity* c_N2R;

    /*------------------------------------------------------------------------*/
    /** \brief Constructor
     *
     * \param[in] ACH the connectivty helper that calls this functor
     * \param[in] AM  the mesh it works on
     */
    BuildR2R_byNFunctor(kmds::GrowingView<kmds::TCellID>* cellIDsAccessor_, ConnectivityHelper* ACH, Mesh* AM, const kmds::Connectivity* Ac_N2R, Kokkos::View<TCellID * [40]> AC2C, Kokkos::View<int*> AC2NB)
            : cellIDsAccessor(cellIDsAccessor_)
            , chelper(ACH)
            , m(AM)
            , c_N2R(Ac_N2R)
            , c2c(AC2C)
            , c2nb(AC2NB)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(TCellID i) const
    {
            TCellID r_id = cellIDsAccessor->get(i);
            Region r = m->getRegion(r_id);
            Kokkos::View<TCellID*> n;
            r.nodeIds(n);

            std::set<TCellID> regionsSet;

            for (auto i_n = 0; i_n < n.size(); i_n++) {
                    TCellID n_id = n(i_n);

                    Kokkos::View<TCellID*> regions_tmp;
                    c_N2R->get(n_id, regions_tmp);

                    for(int i_r=0; i_r<regions_tmp.size(); i_r++) {
                            regionsSet.insert(regions_tmp(i_r));
                    }
            }
            regionsSet.erase(r_id);

            c2nb(r_id)  = regionsSet.size();

            int index = 0;
            for(auto r_tmp: regionsSet) {
                   c2c(r_id, index) = r_tmp;
                   index++;
            }
    }
};
/*------------------------------------------------------------------------*/
//struct BuildR2R_byNFunctor
//{
//    kmds::GrowingView<kmds::TCellID>* cellIDsAccessor;
//    ConnectivityHelper* chelper;
//    Mesh* m;
//    Kokkos::View<TCellID * [40]> c2c;
//    Kokkos::View<int*> c2nb;
//
//    /*------------------------------------------------------------------------*/
//    /** \brief Constructor
//     *
//     * \param[in] ACH the connectivty helper that calls this functor
//     * \param[in] AM  the mesh it works on
//     */
//    BuildR2R_byNFunctor(kmds::GrowingView<kmds::TCellID>* cellIDsAccessor_, ConnectivityHelper* ACH, Mesh* AM, Kokkos::View<TCellID * [40]> AC2C, Kokkos::View<int*> AC2NB)
//            : cellIDsAccessor(cellIDsAccessor_)
//            , chelper(ACH)
//            , m(AM)
//            , c2c(AC2C)
//            , c2nb(AC2NB)
//    {
//    }
//
//    KOKKOS_INLINE_FUNCTION
//    void
//    operator()(int i) const
//    {
//            Region r = m->getRegion(cellIDsAccessor->get(i));
//            Kokkos::View<TCellID*> n;
//            r.nodeIds(n);
//            for (auto i_n = 0; i_n < n.size(); i_n++) {
//                    TCellID n_id = n(i_n);
//                    int top = Kokkos::atomic_fetch_add(&n2nb(n_id), 1);
//                    n2c(n_id, top) = r.id;
//            }
//    }
//};
/*------------------------------------------------------------------------*/
Connectivity::Connectivity(Mesh* AMesh, const EMeshDefinition AT, const TCellID ACapacity)
 : m_mesh(AMesh)
 , m_type(AT)
 , m_first("C_FIRST", ACapacity)
 , m_nb_cells("C_NB", ACapacity)
 , m_C2C("C_TAB", 5 * ACapacity)
 , m_top(0)
 , m_top_C2C(0)
{
}
/*------------------------------------------------------------------------*/
Connectivity::~Connectivity()
{
}
/*------------------------------------------------------------------------*/
EMeshDefinition
Connectivity::getType() const
{
        return m_type;
}
/*------------------------------------------------------------------------*/
void
Connectivity::reinit()
{
        m_top = 0;
        m_top_C2C = 0;
}
/*------------------------------------------------------------------------*/
void
Connectivity::set(const TCellID AKey, const Kokkos::View<TCellID*>& AV)
{
        TSize nb = AV.size();
        TSize pos = Kokkos::atomic_fetch_add(&m_top_C2C, nb);
        m_first(AKey) = pos;
        m_nb_cells(AKey) = nb;
        for (auto i = 0; i < nb; i++) {
                m_C2C(pos + i) = AV(i);
        }
}
/*------------------------------------------------------------------------*/
void
Connectivity::set(const TCellID AKey, const TCellID* AV, const TSize ASize)
{
        TSize nb = ASize;
        TSize pos = Kokkos::atomic_fetch_add(&m_top_C2C, nb);
        m_first(AKey) = pos;
        m_nb_cells(AKey) = nb;
        for (auto i = 0; i < nb; i++) {
                m_C2C(pos + i) = AV[i];
        }
}
/*------------------------------------------------------------------------*/
void
Connectivity::get(const TCellID AKey, Kokkos::View<TCellID*>& AV) const
{
        TSize first = m_first(AKey);
        TSize nb = m_nb_cells(AKey);
        AV = Kokkos::subview(m_C2C, std::make_pair(first, first + nb));
}
/*------------------------------------------------------------------------*/
TSize
Connectivity::getAbsoluteCapacity() const
{
        return m_C2C.extent(0);
}
/*------------------------------------------------------------------------*/
TCellID
Connectivity::getCapacity() const
{
        return m_first.extent(0);
}
/*------------------------------------------------------------------------*/
void
Connectivity::setAbsoluteCapacity(const TSize ASize)
{
        TSize prev_capacity = getAbsoluteCapacity();

        if(prev_capacity < ASize) {

                Kokkos::resize(m_C2C, ASize);

                // We must assign NULLId to new empty items
                TSize from = prev_capacity;
                TSize nb = ASize - prev_capacity;
                Kokkos::View<TCellID *> slice = Kokkos::subview(m_C2C, std::make_pair(from, prev_capacity + nb));
                Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const TSize i) { slice(i) = NullID; });
        }
}
/*------------------------------------------------------------------------*/
void
Connectivity::setCapacity(const TCellID ASize)
{
        TSize prev_capacity = getCapacity();

        if(prev_capacity < ASize) {

                Kokkos::resize(m_first, ASize);
                Kokkos::resize(m_nb_cells, ASize);

                // We must assign NULLId to new empty items
                TSize from = prev_capacity;
                TSize nb = ASize - prev_capacity;
                Kokkos::View<TSize *> slice0 = Kokkos::subview(m_first, std::make_pair(from, prev_capacity + nb));
                Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const TCellID i) { slice0(i) = NullTSize; });
                Kokkos::View<int *> slice1 = Kokkos::subview(m_nb_cells, std::make_pair(from, prev_capacity + nb));
                Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const TCellID i) { slice1(i) = 0; });
        }
}
/*------------------------------------------------------------------------*/
std::string
Connectivity::getName(const EMeshDefinition AD)
{
        std::string s = "undefined";
        switch (AD) {
        case (N2N):
                s = "N2N";
                break;
        case (N2E):
                s = "N2E";
                break;
        case (N2F):
                s = "N2F";
                break;
        case (N2R):
                s = "N2R";
                break;
        case (E2E):
                s = "E2E";
                break;
        case (E2F):
                s = "E2F";
                break;
        case (E2R):
                s = "E2R";
                break;
        case (F2E):
                s = "F2E";
                break;
        case (F2F):
                s = "F2F";
                break;
        case (F2F_byN):
                s = "F2F_byN";
                break;
        case (F2R):
                s = "F2R";
                break;
        case (R2E):
                s = "R2E";
                break;
        case (R2F):
                s = "R2F";
                break;
        case (R2R):
                s = "R2R";
                break;
        case (R2R_byN):
                s = "R2R_byN";
                break;
        default:
                throw KException("Undefined connectivity type");
        };
        return s;
}

void Connectivity::setTop(TSize ASize){m_top = ASize;}
void Connectivity::setTop_C2C(TSize ASize){m_top_C2C = ASize;}
Kokkos::View<TSize*>& Connectivity::getFirst()
{
        return m_first;
}
Kokkos::View<int*>& Connectivity::getNbCells()
{
        return m_nb_cells;
}
Kokkos::View<TCellID*>& Connectivity::getC2C()
{
        return m_C2C;
}
/*------------------------------------------------------------------------*/
ConnectivityHelper::ConnectivityHelper(Mesh* AM)
 : m_mesh(AM)
{
}
/*------------------------------------------------------------------------*/
ConnectivityHelper::~ConnectivityHelper()
{
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildN2E()
{
        TSize nb_nodes = m_mesh->getNbNodes();
        TSize nodes_capacity = m_mesh->getNodeCapacity();
        // For each node, we pre-allocate 10 items
        Kokkos::View<TCellID * [10]> tmp_N2C("tmp_N2C", nodes_capacity);
        // For each node, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", nodes_capacity);
        // init each top to 0
        Kokkos::parallel_for(nodes_capacity, KOKKOS_LAMBDA(const TCellID i) { tmp_nb_cells(i) = 0; });

        // We go through all the cells to fill tmp_N2C
        TSize nb_cells = m_mesh->getNbEdges();
        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getEdgeIDs(&cellIDsAccessor);
        Kokkos::parallel_for(nb_cells, BuildN2EFunctor(&cellIDsAccessor, this, m_mesh, tmp_N2C, tmp_nb_cells));

        // Finally, for each node, we update its N2E connectivity
        Connectivity* cn2c = m_mesh->getConnectivity(N2E);

        // We check if we have enough space to put connectivity in N2C
        TSize nb_n2c = 0;

        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
        m_mesh->getNodeIDs(&nodeIDsAccessor);
        Kokkos::parallel_reduce(nb_nodes, KOKKOS_LAMBDA(const TCellID i, TSize& ls) { ls = ls + tmp_nb_cells(nodeIDsAccessor.get(i)); }, nb_n2c);

        TSize c2nc_size = cn2c->getAbsoluteCapacity();
        if (nb_n2c > c2nc_size) {
                // We need to upgrade the capacity of c2nc
                cn2c->setAbsoluteCapacity(nb_n2c + 10);
        }

        // We reinit the container N2C
        cn2c->reinit();
        Kokkos::parallel_for(nb_nodes, KOKKOS_LAMBDA(const TCellID i) {
            // we udpate the N2C relation for each node
            Node n = m_mesh->getNode(nodeIDsAccessor.get(i));
            // we get the number of cells adjacent to n
            int nb_ad_c = tmp_nb_cells(n.id);
            Kokkos::View<TCellID*> view = Kokkos::subview(tmp_N2C, n.id,std::make_pair(0, nb_ad_c));
            n.setEdges(view);
        });
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildN2F()
{
        TSize nb_nodes = m_mesh->getNbNodes();
        TSize nodes_capacity = m_mesh->getNodeCapacity();
        // For each node, we pre-allocate 20 items
        Kokkos::View<TCellID * [30]> tmp_N2C("tmp_N2C", nodes_capacity);
        // For each node, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", nodes_capacity);
        // init each top to 0
        Kokkos::parallel_for(nodes_capacity, KOKKOS_LAMBDA(const int i) { tmp_nb_cells(i) = 0; });

        // We go through all the cells to fill tmp_N2C
        TSize nb_cells = m_mesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getFaceIDs_dummy(&cellIDsAccessor);
        Kokkos::parallel_for(nb_cells, BuildN2FFunctor(&cellIDsAccessor, this, m_mesh, tmp_N2C, tmp_nb_cells));

        // Finally, for each node, we update its N2F connectivity
        Connectivity* cn2c = m_mesh->getConnectivity(N2F);

        // We check if we have enough space to put connectivity in N2C
        TSize nb_n2c = 0;

        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
        m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);
        Kokkos::parallel_reduce(nb_nodes, KOKKOS_LAMBDA(const TCellID i, TSize& ls) { ls = ls + tmp_nb_cells(nodeIDsAccessor.get(i)); }, nb_n2c);

        TSize c2nc_size = cn2c->getAbsoluteCapacity();
        if (nb_n2c > c2nc_size) {
                // We need to upgrade the capacity of c2nc
                cn2c->setAbsoluteCapacity(nb_n2c + 20);
        }

        // We reinit the container N2C
        cn2c->reinit();
        Kokkos::parallel_for(nb_nodes, KOKKOS_LAMBDA(const TCellID i) {
            // we udpate the N2C relation for each node
            Node n = m_mesh->getNode(nodeIDsAccessor.get(i));
            // we get the number of cells adjacent to n
            int nb_ad_c = tmp_nb_cells(n.id);
            Kokkos::View<TCellID*> view = Kokkos::subview(tmp_N2C, n.id,std::make_pair(0, nb_ad_c));
            n.setFaces(view);
        });
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildN2R()
{
        TSize nb_nodes = m_mesh->getNbNodes();
        TSize nodes_capacity = m_mesh->getNodeCapacity();
        // For each node, we pre-allocate 20 items
        Kokkos::View<TCellID * [20]> tmp_N2C("tmp_N2C", nodes_capacity);
        // For each node, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", nodes_capacity);
        // init each top to 0
        Kokkos::parallel_for(nodes_capacity, KOKKOS_LAMBDA(const TCellID i) { tmp_nb_cells(i) = 0; });

        // We go through all the cells to fill tmp_N2C
        TSize nb_cells = m_mesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getRegionIDs_dummy(&cellIDsAccessor);
        Kokkos::parallel_for(nb_cells, BuildN2RFunctor(&cellIDsAccessor, this, m_mesh, tmp_N2C, tmp_nb_cells));

        // Finally, for each node, we update its N2R connectivity
        Connectivity* cn2c = m_mesh->getConnectivity(N2R);

        // We check if we have enough space to put connectivity in N2C
        TSize nb_n2c = 0;

        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
        m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);
        Kokkos::parallel_reduce(nb_nodes, KOKKOS_LAMBDA(const TCellID i, TSize& ls) { ls = ls + tmp_nb_cells(nodeIDsAccessor.get(i)); }, nb_n2c);

        TSize c2nc_size = cn2c->getAbsoluteCapacity();
        if (nb_n2c > c2nc_size) {
                // We need to upgrade the capacity of c2nc
                cn2c->setAbsoluteCapacity(nb_n2c + 20);
        }

        // We reinit the container N2C
        cn2c->reinit();
        Kokkos::parallel_for(nb_nodes, KOKKOS_LAMBDA(const TCellID i) {
            // we udpate the N2C relation for each node
            Node n = m_mesh->getNode(nodeIDsAccessor.get(i));
            // we get the number of cells adjacent to n
            int nb_ad_c = tmp_nb_cells(n.id);
            Kokkos::View<TCellID*> view = Kokkos::subview(tmp_N2C, n.id,std::make_pair(0, nb_ad_c));
            n.setRegions(view);
//            kmds::TCellID cells[20];
//            for(int ii=0; ii<nb_ad_c; ii++) {
//                    cells[ii] = tmp_N2C(n.id,ii);
//            }
//            n.setRegions(cells, nb_ad_c);

        });
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildN2R_variant_0()
{
    TSize nb_nodes = m_mesh->getNbNodes();
    TSize nodes_capacity = m_mesh->getNodeCapacity();
    // For each node, we pre-allocate 20 items
    Kokkos::View<TCellID * [20]> tmp_N2C("tmp_N2C", nodes_capacity);
    // For each node, we keep in mind how many cells have been added
    Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", nodes_capacity);
    // init each top to 0
    Kokkos::parallel_for(nodes_capacity, KOKKOS_LAMBDA(const int i) { tmp_nb_cells(i) = 0; });

    // We go through all the cells to fill tmp_N2C
    TSize nb_cells = m_mesh->getNbRegions();
//    kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
//    m_mesh->getRegionIDs_dummy(&cellIDsAccessor);
    Kokkos::parallel_for(nb_cells, BuildN2RFunctor_variant_0(this, m_mesh, tmp_N2C, tmp_nb_cells));

    // Finally, for each node, we update its N2R connectivity
    Connectivity* cn2c = m_mesh->getConnectivity(N2R);

    // We check if we have enough space to put connectivity in N2C
    TSize nb_n2c = 0;

//    kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
//    m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);
//    Kokkos::parallel_reduce(nb_nodes, KOKKOS_LAMBDA(const int i, int& ls) { ls = ls + tmp_nb_cells(nodeIDsAccessor.get(i)); }, nb_n2c);
        Kokkos::parallel_reduce(nb_nodes, KOKKOS_LAMBDA(const TCellID i, TSize& ls) { ls = ls + tmp_nb_cells(i); }, nb_n2c);

    TSize c2nc_size = cn2c->getAbsoluteCapacity();
    if (nb_n2c > c2nc_size) {
        // We need to upgrade the capacity of c2nc
        cn2c->setAbsoluteCapacity(nb_n2c + 20);
    }

    // We reinit the container N2C
    cn2c->reinit();
    Kokkos::parallel_for(nb_nodes, KOKKOS_LAMBDA(const TCellID i) {
        // we udpate the N2C relation for each node
//        Node n = m_mesh->getNode(nodeIDsAccessor.get(i));
        Node n = m_mesh->getNode(i);
        // we get the number of cells adjacent to n
        int nb_ad_c = tmp_nb_cells(n.id);
        Kokkos::View<TCellID*> view = Kokkos::subview(tmp_N2C, n.id,std::make_pair(0, nb_ad_c));
        n.setRegions(view);
//            kmds::TCellID cells[20];
//            for(int ii=0; ii<nb_ad_c; ii++) {
//                    cells[ii] = tmp_N2C(n.id,ii);
//            }
//            n.setRegions(cells, nb_ad_c);

    });
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildN2R_variant_1()
{
        Connectivity* cn2c = m_mesh->getConnectivity(N2R);

        // We reinit the container N2C
        cn2c->reinit();

        TSize nb_nodes = m_mesh->getNbNodes();
        // For each node, we keep in mind how many cells have been added

        Kokkos::View<TSize*>& tmp_first = cn2c->getFirst();
        Kokkos::resize(tmp_first, nb_nodes+1);


        TSize nb_cells = m_mesh->getNbRegions();
        Kokkos::parallel_for(nb_cells,
                             KOKKOS_LAMBDA(const TCellID i)
                             {
                                 Region r = m_mesh->getRegion(i);
                                 Kokkos::View<TCellID *> n;
                                 r.nodeIds(n);
                                 for (auto i_n = 0; i_n < n.size(); i_n++) {
                                         TCellID n_id = n(i_n);
                                         TSize top = Kokkos::atomic_fetch_add(&tmp_first(n_id), 1);
                                 }
                             });

        Kokkos::parallel_scan(nb_nodes+1,
                              KOKKOS_LAMBDA(const TCellID& i, kmds::TSize& upd, const bool& final)
                              {
                                  const kmds::TSize val_i = tmp_first(i);
                                  if(final) {
                                          tmp_first(i) = upd;
                                  }

                                  upd += val_i;
                              });

        cn2c->setTop(nb_nodes);
        cn2c->setTop_C2C(tmp_first(nb_nodes));

        Kokkos::View<int*>& tmp_nb_cells = cn2c->getNbCells();
        Kokkos::resize(tmp_nb_cells, nb_nodes+1);

        Kokkos::View<TCellID *>& tmp_c2c = cn2c->getC2C();
        Kokkos::resize(tmp_c2c, tmp_first(nb_nodes));

        Kokkos::parallel_for(nb_cells,
                             KOKKOS_LAMBDA(const TCellID i)
                             {
                                 Region r = m_mesh->getRegion(i);
                                 Kokkos::View<TCellID *> n;
                                 r.nodeIds(n);
                                 for (auto i_n = 0; i_n < n.size(); i_n++) {
                                         TCellID n_id = n(i_n);
                                         TSize top = tmp_first(n_id) + Kokkos::atomic_fetch_add(&tmp_nb_cells(n_id), 1);
                                         tmp_c2c(top) = r.id;
                                 }
                             });

}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildFandF2R()
{
        Kokkos::Timer timer;
        timer.reset();

        TSize nb_cells = m_mesh->getNbRegions();

        // TODO allocate F container
        m_mesh->updateFaceCapacity(nb_cells*8);


        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getRegionIDs_dummy(&cellIDsAccessor);

        struct facecellcell {
                TCellID fid_;
                TCellID cid_[2];
        };

        std::map<FakeFace, facecellcell> ff2id;

        std::cout<<"buildFandF2R preptime "<<timer.seconds()<<std::endl;
        timer.reset();

        for(TCellID i=0; i<nb_cells; i++) {
                TCellID rid = cellIDsAccessor.get(i);
                Region r = m_mesh->getRegion(rid);
                std::vector<FakeFace> fakefaces = r.getFakeFaces();

                for(auto ff: fakefaces) {

                        if(ff2id.find(ff) == ff2id.end()) {
                                // create a new face
                                TCellID fid = m_mesh->newFace(ff);
                                facecellcell fcc;
                                fcc.fid_ = fid;
                                fcc.cid_[0] = rid;
                                fcc.cid_[1] = NullID;
                                ff2id.emplace(ff, fcc);
                        } else {
                                // the face already exists
                                ff2id[ff].cid_[1] = rid;
                        }
                }
        }

        std::cout<<"buildFandF2R fillmap "<<timer.seconds()<<std::endl;
        timer.reset();


        // TODO allocate F2R container
        TSize nb_faces = m_mesh->getNbFaces();
        kmds::Connectivity* c_F2R = m_mesh->getConnectivity(F2R);
        c_F2R->setAbsoluteCapacity(nb_faces*2);
        c_F2R->setCapacity(m_mesh->getFaceCapacity());

        std::cout<<"buildFandF2R prepconnect "<<timer.seconds()<<std::endl;
        timer.reset();

//        for(auto ff: ff2id) {
//                if(ff.second.cid_[1] == NullID) {
//                        Kokkos::View<TCellID *> cc("CELL0CELL1",1);
//                        cc(0) = ff.second.cid_[0];
//                        c_F2R->set(ff.second.fid_, cc);
//                } else {
//                        Kokkos::View<TCellID *> cc("CELL0CELL1",2);
//                        cc(0) = ff.second.cid_[0];
//                        cc(1) = ff.second.cid_[1];
//                        c_F2R->set(ff.second.fid_, cc);
//                }
//
//        }
        for(auto ff: ff2id) {
                TCellID cc[2];
                cc[0] = ff.second.cid_[0];
                cc[1] = ff.second.cid_[1];
                if(ff.second.cid_[1] == NullID) {
                        c_F2R->set(ff.second.fid_, cc, 1);
                } else {
                        c_F2R->set(ff.second.fid_, cc, 2);
                }
        }

        std::cout<<"buildFandF2R fillconnect "<<timer.seconds()<<std::endl;
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildFandF2R_variant_0()
{
//        Kokkos::Timer timer;
//        timer.reset();

        TSize nb_nodes = m_mesh->getNbNodes();
        TSize nodes_capacity = m_mesh->getNodeCapacity();
        // For each node, we pre-allocate 20 items
        Kokkos::View<TCellID * [50][6]> tmp_N2F("tmp_N2F", nodes_capacity);
        // For each node, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", nodes_capacity);
        // init each top to 0
        Kokkos::parallel_for(nodes_capacity, KOKKOS_LAMBDA(const TCellID i) { tmp_nb_cells(i) = 0; });


        TSize nb_cells = m_mesh->getNbRegions();

        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getRegionIDs_dummy(&cellIDsAccessor);

        Kokkos::parallel_for(nb_cells,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 TCellID rid = cellIDsAccessor.get(i);
                                 Region r = m_mesh->getRegion(rid);
                                 std::vector<FakeFace> fakefaces = r.getFakeFaces();

                                 for(auto ff: fakefaces) {
                                         const std::vector<kmds::TCellID> nids = ff.node_ids();
                                         kmds::TCellID firstid = nids[0];
                                         int top = Kokkos::atomic_fetch_add(&tmp_nb_cells(firstid), 1);

                                         if(top >= 50) {
                                             std::cerr<<"TOP greater than 50 "<<firstid<<" "<<top<<" "<<rid<<std::endl;
                                         }

                                         // first store is the size of the face and the cell id
                                         tmp_N2F(firstid,top,0) = nids.size();
                                         tmp_N2F(firstid,top,1) = rid;

                                         for(int ii=1; ii<nids.size(); ii++) {
                                                 tmp_N2F(firstid, top, ii+1) = nids[ii];
                                         }

                                 }
                             });


        m_mesh->updateFaceCapacity(nb_cells*6);

        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
        m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);


        Kokkos::parallel_for(nb_nodes,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 TCellID nid = nodeIDsAccessor.get(i);

                                 // get a fakeface
                                 for(int ff=0; ff<tmp_nb_cells(nid); ff++) {

                                         // check if this face was already created
                                         const int size = tmp_N2F(nid, ff, 0);
                                         if (size == NullID) {
                                                 continue;
                                         } else {
                                                 // create the face
                                                 TCellID nids[size];
                                                 nids[0] = nid;
                                                 for(int ii=1; ii<size; ii++) {
                                                         nids[ii] = tmp_N2F(nid, ff, ii+1);
                                                 }
                                                 TCellID fid = m_mesh->newFace(nids, size);

                                                 // we store the fid of the new face in place of its size
                                                 tmp_N2F(nid, ff, 0) = fid;

                                                 // look for a similar face to get
                                                 bool found = false;
                                                 for(int jj=ff+1; jj<tmp_nb_cells(nid); jj++) {

                                                         bool found_tmp = true;

                                                         const int size_bis = tmp_N2F(nid, jj, 0);

                                                         // that also takes care of NullID
                                                         if (size != size_bis) {
                                                                 found_tmp = false;
                                                                 continue;
                                                         }

                                                         for (int ii = 2; ii <= size_bis; ii++) {
                                                                 if (tmp_N2F(nid, jj, ii) != tmp_N2F(nid, ff, ii)) {
                                                                         found_tmp = false;
                                                                         break;
                                                                 }
                                                         }


                                                         if (found_tmp) {
                                                                 // we repatriate the other cell id
                                                                 // this 2 index is valid because faces always have at least 3 nodes
                                                                 tmp_N2F(nid, ff, 2) = tmp_N2F(nid, jj, 1);

                                                                 // we invalidate this fake face
                                                                 tmp_N2F(nid, jj, 0) = NullID;

                                                                 // if a face was found no need to check the others
                                                                 // as there can be only at most one similar
                                                                 found = true;
                                                                 break;
                                                         }
                                                 }

                                                 if(!found) {
                                                         // we signal that this ff has only one occurence
                                                         tmp_N2F(nid, ff, 2) = NullID;
                                                 }

                                         }

                                 }
                             });

        TSize nb_faces = m_mesh->getNbFaces();
        kmds::Connectivity* c_F2R = m_mesh->getConnectivity(F2R);
        c_F2R->setAbsoluteCapacity(nb_faces*2);
        c_F2R->setCapacity(m_mesh->getFaceCapacity());


        Kokkos::parallel_for(nb_nodes,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 TCellID nid = nodeIDsAccessor.get(i);

                                 // get a fakeface
                                 for (int ff = 0; ff < tmp_nb_cells(nid); ff++) {
                                         // check its validity
                                         if(tmp_N2F(nid, ff, 0) == NullID) {
                                                 continue;
                                         } else {
                                                 TCellID fid = tmp_N2F(nid, ff, 0);
                                                 TCellID cc[2];
                                                 cc[0] = tmp_N2F(nid, ff, 1);
                                                 cc[1] = tmp_N2F(nid, ff, 2);

                                                 if(cc[1] == NullID) {
                                                         c_F2R->set(fid, cc, 1);
                                                 } else {
                                                         c_F2R->set(fid, cc, 2);
                                                 }
                                         }

                                 }
                             });

//        std::cout<<"buildFandF2R fillconnect "<<timer.seconds()<<std::endl;
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildEandE2F()
{
//        Kokkos::Timer timer;
//        timer.reset();

        TSize nb_cells = m_mesh->getNbFaces();

        // TODO allocate E container
        m_mesh->updateEdgeCapacity(nb_cells*5);


        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getFaceIDs_dummy(&cellIDsAccessor);

        struct edgecellcell {
            TCellID eid_;
            TCellID cid_[2];
        };

        std::map<FakeEdge, edgecellcell> fe2id;

//        std::cout<<"buildEandE2F preptime "<<timer.seconds()<<std::endl;
//        timer.reset();

        for(TCellID i=0; i<nb_cells; i++) {
                TCellID fid = cellIDsAccessor.get(i);
                Face f = m_mesh->getFace(fid);
                std::vector<FakeEdge> fakeedges = f.getFakeEdges();

                for(auto fe: fakeedges) {

                        if(fe2id.find(fe) == fe2id.end()) {
                                // create a new edge
                                TCellID eid = m_mesh->newEdge(fe.first(),fe.second());
                                edgecellcell ecc;
                                ecc.eid_ = eid;
                                ecc.cid_[0] = fid;
                                ecc.cid_[1] = NullID;
                                fe2id.emplace(fe, ecc);
                        } else {
                                // the edge already exists
                                fe2id[fe].cid_[1] = fid;
                        }
                }
        }
//        std::cout<<"buildEandE2F filledmap "<<timer.seconds()<<std::endl;
//        timer.reset();

        // TODO allocate E2F container
        TSize nb_edges = m_mesh->getNbEdges();
        kmds::Connectivity* c_E2F = m_mesh->getConnectivity(E2F);
        c_E2F->setAbsoluteCapacity(nb_edges*2);
        c_E2F->setCapacity(m_mesh->getEdgeCapacity());

//        std::cout<<"buildEandE2F prepconnect "<<timer.seconds()<<std::endl;
//        timer.reset();

//        for(auto fe: fe2id) {
//                if(fe.second.cid_[1] == NullID) {
//                        Kokkos::View<TCellID *> cc("CELL0CELL1",1);
//                        cc(0) = fe.second.cid_[0];
//                        c_E2F->set(fe.second.eid_, cc);
//                } else {
//                        Kokkos::View<TCellID *> cc("CELL0CELL1",2);
//                        cc(0) = fe.second.cid_[0];
//                        cc(1) = fe.second.cid_[1];
//                        c_E2F->set(fe.second.eid_, cc);
//                }
//
//        }
        for(auto fe: fe2id) {
                TCellID cc[2];
                cc[0] = fe.second.cid_[0];
                cc[1] = fe.second.cid_[1];
                if(fe.second.cid_[1] == NullID) {
                        c_E2F->set(fe.second.eid_, cc, 1);
                } else {
                        c_E2F->set(fe.second.eid_, cc, 2);
                }
        }
//        std::cout<<"buildEandE2F filled connect "<<timer.seconds()<<std::endl;
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildEandE2F_2D_variant_0()
{
//        Kokkos::Timer timer;
//        timer.reset();

        TSize nb_nodes = m_mesh->getNbNodes();
        TSize nodes_capacity = m_mesh->getNodeCapacity();
        // For each node, we pre-allocate 20 items
        Kokkos::View<TCellID * [30][3]> tmp_N2E("tmp_N2E", nodes_capacity);
        // For each node, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", nodes_capacity);
        // init each top to 0
        Kokkos::parallel_for(nodes_capacity, KOKKOS_LAMBDA(const TCellID i) { tmp_nb_cells(i) = 0; });


        TSize nb_cells = m_mesh->getNbFaces();

        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getFaceIDs_dummy(&cellIDsAccessor);

        Kokkos::parallel_for(nb_cells,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 TCellID fid = cellIDsAccessor.get(i);
                                 Face f = m_mesh->getFace(fid);
                                 std::vector<FakeEdge> fakeedges = f.getFakeEdges();

                                 for(auto fe: fakeedges) {

                                         const TCellID firstid = fe.first();
                                         const TCellID lastid = fe.second();

                                         const int top = Kokkos::atomic_fetch_add(&tmp_nb_cells(firstid), 1);

                                         // first store is the size of the face and the cell id
                                         tmp_N2E(firstid,top,1) = fid;
                                         tmp_N2E(firstid,top,2) = lastid;
                                 }
                             });


        m_mesh->updateEdgeCapacity(nb_cells*4);

        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
        m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);


        Kokkos::parallel_for(nb_nodes,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 TCellID nid = nodeIDsAccessor.get(i);

                                 // get a fakeedge
                                 for(int fe=0; fe<tmp_nb_cells(nid); fe++) {

                                         // check if this edge was already created
                                         if (tmp_N2E(nid, fe, 0) == NullID) {
                                                 continue;
                                         } else {
                                                 // create the edge
                                                 TCellID lastid = tmp_N2E(nid, fe, 2);
                                                 TCellID eid = m_mesh->newEdge(nid, lastid);

                                                 // we store the eid of the new edge
                                                 tmp_N2E(nid, fe, 0) = eid;

                                                 // look for a similar edge to get
                                                 bool found = false;
                                                 for(int jj=fe+1; jj<tmp_nb_cells(nid); jj++) {

                                                         if(tmp_N2E(nid, jj, 2) == lastid) {

                                                                 // we repatriate the other cell id
                                                                 tmp_N2E(nid, fe, 2) = tmp_N2E(nid, jj, 1);

                                                                 // we invalidate this fake edge
                                                                 tmp_N2E(nid, jj, 0) = NullID;

                                                                 found = true;
                                                                 break;
                                                         }
                                                 }

                                                 if(!found) {
                                                         // we signal that this fe has only one occurence
                                                         tmp_N2E(nid, fe, 2) = NullID;
                                                 }

                                         }

                                 }
                             });

        TSize nb_edges = m_mesh->getNbEdges();
        kmds::Connectivity* c_E2F = m_mesh->getConnectivity(E2F);
        c_E2F->setAbsoluteCapacity(nb_edges*2);
        c_E2F->setCapacity(m_mesh->getEdgeCapacity());


        Kokkos::parallel_for(nb_nodes,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 TCellID nid = nodeIDsAccessor.get(i);

                                 // get a fakeedge
                                 for (int fe = 0; fe < tmp_nb_cells(nid); fe++) {
                                         // check its validity
                                         TCellID eid = tmp_N2E(nid, fe, 0);
                                         if(eid == NullID) {
                                                 continue;
                                         } else {
                                                 TCellID cc[2];
                                                 cc[0] = tmp_N2E(nid, fe, 1);
                                                 cc[1] = tmp_N2E(nid, fe, 2);

                                                 if(cc[1] == NullID) {
                                                         c_E2F->set(eid, cc, 1);
                                                 } else {
                                                         c_E2F->set(eid, cc, 2);
                                                 }
                                         }

                                 }
                             });

//        std::cout<<"buildFandF2R fillconnect "<<timer.seconds()<<std::endl;
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildEandE2F_2D_variant_1()
{
//        Kokkos::Timer timer;
//        timer.reset();

        TSize nb_cells = m_mesh->getNbFaces();
        TSize nb_nodes = m_mesh->getNbNodes();

        // max theoretical number of egdes
        TSize maxNbEdgesTh = nb_cells*4;

        m_mesh->updateEdgeCapacity(maxNbEdgesTh);

        kmds::Connectivity* c_E2F = m_mesh->getConnectivity(E2F);
        c_E2F->setAbsoluteCapacity(maxNbEdgesTh*2);
        c_E2F->setCapacity(m_mesh->getEdgeCapacity());

        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getFaceIDs_dummy(&cellIDsAccessor);


        //Kokkos::UnorderedMap<uint64_t, TCellID> e2fid_single(maxNbEdgesTh);
        Kokkos::UnorderedMap<uint64_t, void> e2fid_single(maxNbEdgesTh);
        Kokkos::UnorderedMap<uint64_t, TCellID> e2fid_double(maxNbEdgesTh);


        for(TCellID i=0; i<nb_cells; i++) {
                TCellID fid = cellIDsAccessor.get(i);
                Face f = m_mesh->getFace(fid);
                std::vector<FakeEdge> fakeedges = f.getFakeEdges();

                for(auto fe: fakeedges) {

                        // unique id for the edge, defined as its expression in base maxNbEdgesTh
                        uint64_t eid_unique = fe.second() + fe.first() * nb_nodes;

                        //Kokkos::UnorderedMapInsertResult res = e2fid_single.insert(eid_unique, fid);
                        Kokkos::UnorderedMapInsertResult res = e2fid_single.insert(eid_unique);
//                        bool success = res.success();
                        bool exist = res.existing();
//                        bool fail = res.failed();

                        if(exist) {
                                Kokkos::UnorderedMapInsertResult res_bis = e2fid_double.insert(eid_unique, fid);
                        } else {
                                TCellID eid = m_mesh->newEdge(fe.first(), fe.second());
                                TCellID cc[2];
                                cc[0] = fid;
                                cc[1] = NullID;
                                c_E2F->set(eid, cc, 2);
                        }
                }
        }
//        std::cout<<"buildEandE2F filledmap "<<timer.seconds()<<std::endl;
//        timer.reset();

        TSize nb_edges = m_mesh->getNbEdges();
        kmds::GrowingView<kmds::TCellID> edgeIDsAccessor("EDGEIDSACCESSOR", nb_edges);
        m_mesh->getEdgeIDs_dummy(&edgeIDsAccessor);

        Kokkos::View<TSize*>& tmp_first = c_E2F->getFirst();
        Kokkos::View<int*>& tmp_nb_cells = c_E2F->getNbCells();
        Kokkos::View<TCellID *>& tmp_c2c = c_E2F->getC2C();

        Kokkos::parallel_for(nb_edges,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 Edge e = m_mesh->getEdge(edgeIDsAccessor.get(i));
                                 TCellID ids[2];
                                 e.nodeIds(ids);

                                 uint64_t eid_unique = ids[1] + ids[0] * nb_nodes;

                                 uint32_t index = e2fid_double.find(eid_unique);

                                 if(e2fid_double.valid_at(index)) {
                                         TCellID fid = e2fid_double.value_at(index);
                                         tmp_c2c(tmp_first(e.id) + 1) = fid;
                                 } else {
                                         tmp_nb_cells(e.id) = 1;
                                 }
                             });

//        std::cout<<"buildEandE2F filled connect "<<timer.seconds()<<std::endl;
}

/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildEandE2F_3D()
{
//        Kokkos::Timer timer;
//        timer.reset();

    TSize nb_cells = m_mesh->getNbFaces();
    TSize nbNodes = m_mesh->getNbNodes();

    kmds::GrowingView<kmds::TCellID> faceIDsAccessor("FACEIDSACCESSOR", nb_cells);
    m_mesh->getFaceIDs(&faceIDsAccessor);

    // max theoretical number of egdes
    TSize maxNbEdgesTh = nb_cells*5;

    Kokkos::UnorderedMap<uint64_t, void> fe_set(maxNbEdgesTh);

    // build a fakeedges set
    Kokkos::parallel_for(nb_cells,
                         KOKKOS_LAMBDA(const TCellID i) {
                             Face f = m_mesh->getFace(faceIDsAccessor.get(i));
                             std::vector<FakeEdge> fakeedges = f.getFakeEdges();

                             for (auto fe: fakeedges) {
                                 uint64_t eid_unique = fe.second() + fe.first() * nbNodes;

                                 Kokkos::UnorderedMapInsertResult res = fe_set.insert(eid_unique);
                             }
                         }
    );

    // edges creation and fakeedge to edgeID mapping
    TSize nbEdges = fe_set.size();
    m_mesh->updateEdgeCapacity(nbEdges);

    Kokkos::UnorderedMap<uint64_t, TCellID> fe2eid(nbEdges);

    Kokkos::parallel_for(fe_set.capacity(),
                         KOKKOS_LAMBDA(const TCellID i) {
                             if(fe_set.valid_at(i)) {
                                 uint64_t eid_unique = fe_set.key_at(i);

                                 std::ldiv_t d = std::ldiv(eid_unique, nbNodes);

                                 TCellID nid0 = d.rem;
                                 TCellID nid1 = d.quot;

                                 TCellID eid = m_mesh->newEdge(nid0, nid1);
                                 fe2eid.insert(eid_unique, eid);
                             }
                         }
    );

    Connectivity* c_e2c = m_mesh->getConnectivity(E2F);

    // We reinit the container N2C
    c_e2c->reinit();

    Kokkos::View<TSize*>& tmp_first = c_e2c->getFirst();
    Kokkos::resize(tmp_first, nbEdges+1);

    Kokkos::parallel_for(nb_cells,
                         KOKKOS_LAMBDA(const TCellID i) {
                             Face f = m_mesh->getFace(faceIDsAccessor.get(i));
                             std::vector<FakeEdge> fakeedges = f.getFakeEdges();

                             for (auto fe: fakeedges) {
                                 uint64_t eid_unique = fe.second() + fe.first() * nbNodes;

                                 TSize index = fe2eid.find(eid_unique);
                                 TCellID eid = fe2eid.value_at(index);

                                 TSize top = Kokkos::atomic_fetch_add(&tmp_first(eid), 1);
                             }
                         });

    Kokkos::parallel_scan(nbEdges+1,
                          KOKKOS_LAMBDA(const TCellID& i, TSize& upd, const bool& final)
                          {
                              const kmds::TSize val_i = tmp_first(i);
                              if(final) {
                                  tmp_first(i) = upd;
                              }

                              upd += val_i;
                          });

    c_e2c->setTop(nbEdges);
    c_e2c->setTop_C2C(tmp_first(nbEdges));

    Kokkos::View<int*>& tmp_nb_cells = c_e2c->getNbCells();
    Kokkos::resize(tmp_nb_cells, nbEdges+1);

    Kokkos::View<TCellID *>& tmp_c2c = c_e2c->getC2C();
    Kokkos::resize(tmp_c2c, tmp_first(nbEdges));

    Kokkos::parallel_for(nb_cells,
                         KOKKOS_LAMBDA(const TCellID i)
                         {
                             Face f = m_mesh->getFace(i);
                             std::vector<FakeEdge> fakeedges = f.getFakeEdges();

                             for (auto fe: fakeedges) {
                                 uint64_t eid_unique = fe.second() + fe.first() * nbNodes;

                                 TSize index = fe2eid.find(eid_unique);
                                 TCellID eid = fe2eid.value_at(index);
                                 TSize top = tmp_first(eid) + Kokkos::atomic_fetch_add(&tmp_nb_cells(eid), 1);
                                 tmp_c2c(top) = f.id;
                             }
                         });

}

/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildEandE2R()
{
//        Kokkos::Timer timer;
//        timer.reset();

        TSize nb_cells = m_mesh->getNbRegions();
        TSize nbNodes = m_mesh->getNbNodes();

        kmds::GrowingView<kmds::TCellID> regionIDsAccessor("REGIONIDSACCESSOR", nb_cells);
        m_mesh->getRegionIDs(&regionIDsAccessor);

        // max theoretical number of egdes
        TSize maxNbEdgesTh = nb_cells*12;

        Kokkos::UnorderedMap<uint64_t, void> fe_set(maxNbEdgesTh);

        // build a fakeedges set
        Kokkos::parallel_for(nb_cells,
                             KOKKOS_LAMBDA(const TCellID i) {
                                 Region r = m_mesh->getRegion(regionIDsAccessor.get(i));
                                 std::vector<FakeEdge> fakeedges = r.getFakeEdges();

                                 for (auto fe: fakeedges) {
                                         uint64_t eid_unique = fe.second() + fe.first() * nbNodes;

                                         Kokkos::UnorderedMapInsertResult res = fe_set.insert(eid_unique);
                                 }
                             }
        );

        // edges creation and fakeedge to edgeID mapping
        TSize nbEdges = fe_set.size();
        m_mesh->updateEdgeCapacity(nbEdges);

        Kokkos::UnorderedMap<uint64_t, TCellID> fe2eid(nbEdges);

        Kokkos::parallel_for(fe_set.capacity(),
                             KOKKOS_LAMBDA(const TSize i) {
                                 if(fe_set.valid_at(i)) {
                                         uint64_t eid_unique = fe_set.key_at(i);

                                         std::ldiv_t d = std::ldiv(eid_unique, nbNodes);

                                         TCellID nid0 = d.rem;
                                         TCellID nid1 = d.quot;

                                         TCellID eid = m_mesh->newEdge(nid0, nid1);
                                         fe2eid.insert(eid_unique, eid);
                                 }
                             }
        );

        Connectivity* c_e2c = m_mesh->getConnectivity(E2R);

        // We reinit the container N2C
        c_e2c->reinit();

        Kokkos::View<TSize*>& tmp_first = c_e2c->getFirst();
        Kokkos::resize(tmp_first, nbEdges+1);

        Kokkos::parallel_for(nb_cells,
                             KOKKOS_LAMBDA(const TCellID i) {
                Region r = m_mesh->getRegion(regionIDsAccessor.get(i));
                std::vector<FakeEdge> fakeedges = r.getFakeEdges();

                for (auto fe: fakeedges) {
                        uint64_t eid_unique = fe.second() + fe.first() * nbNodes;

                        TSize index = fe2eid.find(eid_unique);
                        TCellID eid = fe2eid.value_at(index);

                        TSize top = Kokkos::atomic_fetch_add(&tmp_first(eid), 1);
                }
        });

        Kokkos::parallel_scan(nbEdges+1,
                              KOKKOS_LAMBDA(const TCellID& i, TSize& upd, const bool& final)
        {
                const kmds::TSize val_i = tmp_first(i);
                if(final) {
                        tmp_first(i) = upd;
                }

                upd += val_i;
        });

        c_e2c->setTop(nbEdges);
        c_e2c->setTop_C2C(tmp_first(nbEdges));

        Kokkos::View<int*>& tmp_nb_cells = c_e2c->getNbCells();
        Kokkos::resize(tmp_nb_cells, nbEdges+1);

        Kokkos::View<TCellID *>& tmp_c2c = c_e2c->getC2C();
        Kokkos::resize(tmp_c2c, tmp_first(nbEdges));

        Kokkos::parallel_for(nb_cells,
                             KOKKOS_LAMBDA(const TCellID i)
        {
                Region r = m_mesh->getRegion(i);
                std::vector<FakeEdge> fakeedges = r.getFakeEdges();

                for (auto fe: fakeedges) {
                        uint64_t eid_unique = fe.second() + fe.first() * nbNodes;

                        TSize index = fe2eid.find(eid_unique);
                        TCellID eid = fe2eid.value_at(index);
                        TSize top = tmp_first(eid) + Kokkos::atomic_fetch_add(&tmp_nb_cells(eid), 1);
                        tmp_c2c(top) = r.id;
                }
        });

}

/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildF2F_byN(const kmds::Connectivity* c_N2F)
{
        TSize nb_cells = m_mesh->getNbFaces();
        TSize cells_capacity = m_mesh->getFaceCapacity();
        // For each item, we pre-allocate 40 items
        Kokkos::View<TCellID * [20]> tmp_C2C("tmp_C2C", cells_capacity);
        // For each item, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", cells_capacity);
        // init each top to 0
        Kokkos::parallel_for(cells_capacity, KOKKOS_LAMBDA(const TCellID i) { tmp_nb_cells(i) = 0; });

        // We go through all the cells to fill tmp_C2C
        //TSize nb_cells = m_mesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getFaceIDs_dummy(&cellIDsAccessor);
        Kokkos::parallel_for(nb_cells, BuildF2F_byNFunctor(&cellIDsAccessor, this, m_mesh, c_N2F, tmp_C2C, tmp_nb_cells));

        // Finally, for each cell, we update its C2C_byN connectivity
        Connectivity* c_c2c = m_mesh->getConnectivity(F2F_byN);

        // We check if we have enough space to put connectivity in C2C
        TSize nb_c2c = 0;

//        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
//        m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);
        Kokkos::parallel_reduce(nb_cells, KOKKOS_LAMBDA(const TCellID i, TSize& ls) { ls = ls + tmp_nb_cells(cellIDsAccessor.get(i)); }, nb_c2c);

        TSize c2c_size = c_c2c->getAbsoluteCapacity();
        if (nb_c2c > c2c_size) {
                // We need to upgrade the capacity of c2nc
                c_c2c->setAbsoluteCapacity(nb_c2c + 20);
        }

        // We reinit the container C2C
        c_c2c->reinit();
        Kokkos::parallel_for(nb_cells, KOKKOS_LAMBDA(const TCellID i) {
            // we udpate the C2C relation for each cell
            TCellID cid = cellIDsAccessor.get(i);
            // we get the number of cells adjacent to c
            int nb_ad_c = tmp_nb_cells(cid);
            Kokkos::View<TCellID*> view = Kokkos::subview(tmp_C2C, cid,std::make_pair(0, nb_ad_c));
            c_c2c->set(cid, view);
        });
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildR2R_byN(const kmds::Connectivity* c_N2R)
{
        TSize nb_cells = m_mesh->getNbRegions();
        TSize cells_capacity = m_mesh->getRegionCapacity();
        // For each item, we pre-allocate 40 items
        Kokkos::View<TCellID * [40]> tmp_C2C("tmp_C2C", cells_capacity);
        // For each item, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", cells_capacity);
        // init each top to 0
        Kokkos::parallel_for(cells_capacity, KOKKOS_LAMBDA(const TCellID i) { tmp_nb_cells(i) = 0; });

        // We go through all the cells to fill tmp_C2C
        //TSize nb_cells = m_mesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
        m_mesh->getRegionIDs_dummy(&cellIDsAccessor);
        Kokkos::parallel_for(nb_cells, BuildR2R_byNFunctor(&cellIDsAccessor, this, m_mesh, c_N2R, tmp_C2C, tmp_nb_cells));

        // Finally, for each cell, we update its C2C_byN connectivity
        Connectivity* c_c2c = m_mesh->getConnectivity(R2R_byN);

        // We check if we have enough space to put connectivity in C2C
        TSize nb_c2c = 0;

        //kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
        //m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);
        Kokkos::parallel_reduce(nb_cells, KOKKOS_LAMBDA(const TCellID i, TSize & ls) { ls = ls + tmp_nb_cells(cellIDsAccessor.get(i)); }, nb_c2c);

        TSize c2c_size = c_c2c->getAbsoluteCapacity();
        if (nb_c2c > c2c_size) {
                // We need to upgrade the capacity of c2nc
                c_c2c->setAbsoluteCapacity(nb_c2c + 40);
        }

        // We reinit the container C2C
        c_c2c->reinit();
        Kokkos::parallel_for(nb_cells, KOKKOS_LAMBDA(const TCellID i) {
            // we udpate the C2C relation for each cell
            TCellID cid = cellIDsAccessor.get(i);
            // we get the number of cells adjacent to c
            int nb_ad_c = tmp_nb_cells(cid);
            Kokkos::View<TCellID*> view = Kokkos::subview(tmp_C2C, cid,std::make_pair(0, nb_ad_c));
            c_c2c->set(cid, view);
        });
}
/*------------------------------------------------------------------------*/
void
ConnectivityHelper::buildN2N(const EMeshDefinition ACellType)
{
        std::set<FakeEdge> allFakeEdges;

        switch(ACellType) {

                case F :
                {
                        TSize nb_cells = m_mesh->getNbFaces();
                        kmds::GrowingView<TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
                        m_mesh->getFaceIDs_dummy(&cellIDsAccessor);

                        for(int i=0; i<nb_cells; i++) {
                                Face f = m_mesh->getFace(cellIDsAccessor.get(i));
                                std::vector<FakeEdge> fakeEdges = f.getFakeEdges();

                                for(auto fe: fakeEdges) {
                                        allFakeEdges.insert(fe);
                                }
                        }
                }
                break;
                case R :
                {
                        TSize nb_cells = m_mesh->getNbRegions();
                        kmds::GrowingView<TCellID> cellIDsAccessor("CELLIDSACCESSOR", nb_cells);
                        m_mesh->getRegionIDs_dummy(&cellIDsAccessor);

                        for(int i=0; i<nb_cells; i++) {
                                Region r = m_mesh->getRegion(cellIDsAccessor.get(i));
                                std::vector<FakeEdge> fakeEdges = r.getFakeEdges();

                                for(auto fe: fakeEdges) {
                                        allFakeEdges.insert(fe);
                                }
                        }
                }
                break;
                default:
                        throw KException("ConnectivityHelper::buildN2N cell type not recognised.");
        }

        TSize nb_nodes = m_mesh->getNbNodes();
        TSize nodes_capacity = m_mesh->getNodeCapacity();

        // For each node, we pre-allocate 20 items
        Kokkos::View<TCellID * [20]> tmp_N2C("tmp_N2C", nodes_capacity);
        // For each node, we keep in mind how many cells have been added
        Kokkos::View<int*> tmp_nb_cells("tmp_nb_cells", nodes_capacity);
        // init each top to 0
        Kokkos::parallel_for(nodes_capacity, KOKKOS_LAMBDA(const TCellID i) { tmp_nb_cells(i) = 0; });

        for(auto fe: allFakeEdges) {
                TCellID id0 = fe.first();
                TCellID id1 = fe.second();

                tmp_N2C(id0, tmp_nb_cells(id0)) = id1;
                tmp_N2C(id1, tmp_nb_cells(id1)) = id0;


                tmp_nb_cells(id0) = tmp_nb_cells(id0) + 1;
                tmp_nb_cells(id1) = tmp_nb_cells(id1) + 1;
        }

        // Finally, for each node, we update its N2R connectivity
        Connectivity* cn2c = m_mesh->getConnectivity(N2N);

        // We check if we have enough space to put connectivity in N2C
        TSize nb_n2c = 0;

        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODEIDSACCESSOR", nb_nodes);
        m_mesh->getNodeIDs_dummy(&nodeIDsAccessor);
        Kokkos::parallel_reduce(nb_nodes, KOKKOS_LAMBDA(const TCellID i, TSize& ls) { ls = ls + tmp_nb_cells(nodeIDsAccessor.get(i)); }, nb_n2c);

        TSize c2nc_size = cn2c->getAbsoluteCapacity();
        if (nb_n2c > c2nc_size) {
                // We need to upgrade the capacity of c2nc
                cn2c->setAbsoluteCapacity(nb_n2c + 20);
        }

        // We reinit the container N2C
        cn2c->reinit();
        Kokkos::parallel_for(nb_nodes, KOKKOS_LAMBDA(const TCellID i) {
            // we udpate the N2C relation for each node
            TCellID nid = nodeIDsAccessor.get(i);
            // we get the number of cells adjacent to n
            int nb_ad_c = tmp_nb_cells(nid);
            Kokkos::View<TCellID*> view = Kokkos::subview(tmp_N2C, nid,std::make_pair(0, nb_ad_c));
            cn2c->set(nid, view);
        });
}
/*------------------------------------------------------------------------*/
