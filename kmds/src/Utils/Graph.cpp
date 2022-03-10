/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * Graph.cpp
 *
 *  Created on: 6 mars 2018
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include "KM/Utils/Graph.h"
/*----------------------------------------------------------------------------*/
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
#include <bitset>
#include <random>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
    TCellID
    Graph::getNbVec() const
    {
        return m_nbvec;
    }
/*----------------------------------------------------------------------------*/
    TCellID
    Graph::getNbEdges() const
    {
        int s = 0;
        Kokkos::parallel_reduce(this->getNbVec(),
                                KOKKOS_LAMBDA(const int i, int& ls) {
                                    Kokkos::View<kmds::TCellID *> neighbors;
                                    int nb = this->getNeighbors(i, neighbors);
                                    ls = ls + nb;
                                },
                                s);
        return s;
    }
/*----------------------------------------------------------------------------*/
    void
    Graph::setNeighbors(const TCellID AId, const std::set<kmds::TCellID>& AN)
    {
        m_nbneighbors(AId) = AN.size();
        int index=0;
        for(auto n: AN) {
            m_neighbors(AId,index) = n;
            index++;
        }
    }
/*----------------------------------------------------------------------------*/
    void
    Graph::setNeighbors(const TCellID AId, const std::vector<kmds::TCellID>& AN)
    {
        m_nbneighbors(AId) = AN.size();
        int index=0;
        for(auto n: AN) {
            m_neighbors(AId,index) = n;
            index++;
        }
    }
/*----------------------------------------------------------------------------*/
    TCellID
    Graph::getNeighbors(const TCellID AId, Kokkos::View<TCellID*>& AN) const
    {
        int nb = m_nbneighbors(AId);
        AN = Kokkos::subview(m_neighbors, AId, std::make_pair(0, nb));

        return nb;
    }
    /*----------------------------------------------------------------------------*/
    bool
    Graph::isNeighbor(const TCellID AId0, const TCellID AId1) const
    {
        int nb = m_nbneighbors(AId0);
        Kokkos::View<TCellID*> AN;
        AN = Kokkos::subview(m_neighbors, AId0, std::make_pair(0, nb));

        for(int i_n=0; i_n<nb; i_n++) {
            if(AId1 == AN[i_n]) {
                return true;
            }
        }

        return false;
    }
/*----------------------------------------------------------------------------*/
    void
    Graph::getIndependentSet(kmds::GrowingView<kmds::TCellID> *ASelection) const
    {
        // assigned vertices colors, 0 is the unassigned color

        // THREAD atomic trait on this view
        Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic> > coloring("COLORING", this->getNbVec());


        struct ColorCell
        {
            const Graph* gr;
            Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic> >* coloring;
            kmds::GrowingView<kmds::TCellID> *selection;

            ColorCell(const Graph* gr_, Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic> >* coloring_, kmds::GrowingView<kmds::TCellID>* selection_)
                    : gr(gr_)
                    , coloring(coloring_)
                    , selection(selection_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                TCellID id = selection->get(i);

                Kokkos::View<kmds::TCellID *> neighbors;
                int nb = gr->getNeighbors(id, neighbors);

                // check the colors already assigned to neighbors
                std::bitset<MaxNbColors+1> forbidden(false);

                for(int ineighbor=0; ineighbor<nb; ineighbor++) {

                    // THREAD atomic read
                    int color = (*coloring)[neighbors[ineighbor]];

                    forbidden[color] = true;

                }

                // select an allowed color; we begin from 1 because 0 counts as unassigned
                int freeCol;
                for(int icol=1; icol<=Graph::MaxNbColors; icol++) {
                    if(!forbidden[icol]) {

                        // THREAD atomic write
                        (*coloring)[id] = icol;
                        break;
                    }
                }
            }
        };

        // first assignement that  can contain conflicts
        kmds::GrowingView<kmds::TCellID> conflictvertices("CONFLICT_VERTICES", this->getNbVec());

        struct selectionPopulate
        {
            kmds::GrowingView<kmds::TCellID> *selection;

            selectionPopulate(kmds::GrowingView<kmds::TCellID>* selection_)
            : selection(selection_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                selection->push_back(i);
            }
        };

        Kokkos::parallel_for(this->getNbVec(), selectionPopulate(&conflictvertices));
        Kokkos::parallel_for(this->getNbVec(), ColorCell(this, &coloring, &conflictvertices));


        // get the conflicts
        struct ConflictSelection
        {
            const Graph *gr;
            Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic> > *coloring;
            kmds::GrowingView<kmds::TCellID> *selection;

            ConflictSelection(const Graph *gr_, Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic> > *coloring_, kmds::GrowingView<kmds::TCellID>* selection_)
                    : gr(gr_)
                    , coloring(coloring_)
                    , selection(selection_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
//                TCellID id = selection->get(i);
                TCellID id = i;

                Kokkos::View<kmds::TCellID *> neighbors;
                int nb = gr->getNeighbors(id, neighbors);

                int ownColor = (*coloring)[id];

                for(int ineighbor=0; ineighbor<nb; ineighbor++) {

                    TCellID neighborID = neighbors[ineighbor];

                    if((*coloring)[neighborID] == ownColor) {

                        // we keep only one of the two, the one with the highest id
                        if (neighborID < id) {
                            selection->push_back(id);
                            (*coloring)[id] = 0;
                            break;
                        }
                    }
                }
            }
        };

        // detect and solve the conflicts
        conflictvertices.clear();

        Kokkos::parallel_for(this->getNbVec(), ConflictSelection(this, &coloring, &conflictvertices));

        int nbConflicts = conflictvertices.getNbElems();

        while (nbConflicts > 0) {
            std::cout<<"nbConflicts "<<nbConflicts<<std::endl;

            Kokkos::parallel_for(nbConflicts, ColorCell(this, &coloring, &conflictvertices));

            conflictvertices.clear();

            Kokkos::parallel_for(this->getNbVec(), ConflictSelection(this, &coloring, &conflictvertices));
            nbConflicts =  conflictvertices.getNbElems();

            //exit(0);
        }

        // select vertices assigned color == 1
        Kokkos::parallel_for(this->getNbVec(), KOKKOS_LAMBDA(const int i)
            {
                if(coloring[i] == 1) {
                    ASelection->push_back(i);
                }
            }
        );



    }
/*----------------------------------------------------------------------------*/
    bool
    Graph::checkIndependentSet(kmds::GrowingView<kmds::TCellID>* ASelection) const
    {
        Kokkos::UnorderedMap<kmds::TCellID, void> kmap(this->getNbVec());

        int nbFailed = 0;

        Kokkos::parallel_reduce(ASelection->getNbElems(), KOKKOS_LAMBDA(const int i, int &sum) {

                                    Kokkos::UnorderedMapInsertResult res = kmap.insert(ASelection->get(i));
                                    bool success = res.success();
                                    bool exist = res.existing();
                                    bool fail = res.failed();

                                    if (fail) {
                                        sum++;
                                    }
                                },
                                nbFailed
        );

        int found = 0;
        Kokkos::parallel_reduce(ASelection->getNbElems(), KOKKOS_LAMBDA(const int i, int &res) {

                                    TCellID id = ASelection->get(i);

                                    Kokkos::View<kmds::TCellID *> neighbors;
                                    int nb = this->getNeighbors(id, neighbors);

                                    for(int ineighbor=0; ineighbor<nb; ineighbor++) {
                                        TCellID neighborID = neighbors[ineighbor];

                                        int index = kmap.find(neighborID);
                                        if(kmap.valid_at(index)) {
                                            std::cout<<"found "<<id<<" "<<neighborID<<std::endl;
                                            res++;
                                        }
                                    }

                                },
                                found
        );

        return (found == 0);
    }
/*----------------------------------------------------------------------------*/
    void
    Graph::buildColoring(const int nbMaxColors)
    {
        for(int icolor=0; icolor<nbMaxColors; icolor++) {
            m_colorsets.push_back(kmds::GrowingView<kmds::TCellID> ("COLOR_SET", this->getNbVec()));
        }

        // assigned vertices colors, 0 is the unassigned color

        // THREAD atomic trait on this view
        Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic> > coloring("COLORING", this->getNbVec());


        struct ColorCell
        {
            const Graph* gr;
            Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic> >* coloring;
            kmds::GrowingView<kmds::TCellID> *selection;

            ColorCell(const Graph* gr_, Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic> >* coloring_, kmds::GrowingView<kmds::TCellID>* selection_)
                    : gr(gr_)
                    , coloring(coloring_)
                    , selection(selection_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                TCellID id = selection->get(i);

                Kokkos::View<kmds::TCellID *> neighbors;
                int nb = gr->getNeighbors(id, neighbors);

                // check the colors already assigned to neighbors
                std::bitset<MaxNbColors+1> forbidden(false);

                for(int ineighbor=0; ineighbor<nb; ineighbor++) {

                    // THREAD atomic read
                    int color = (*coloring)[neighbors[ineighbor]];

                    forbidden[color] = true;

                }

                // select an allowed color; we begin from 1 because 0 counts as unassigned
                int freeCol;
                for(int icol=1; icol<=Graph::MaxNbColors; icol++) {
                    if(!forbidden[icol]) {

                        // THREAD atomic write
                        (*coloring)[id] = icol;
                        break;
                    }
                }
            }
        };

        // first assignement that  can contain conflicts
        kmds::GrowingView<kmds::TCellID> conflictvertices("CONFLICT_VERTICES", this->getNbVec());

        struct selectionPopulate
        {
            kmds::GrowingView<kmds::TCellID> *selection;

            selectionPopulate(kmds::GrowingView<kmds::TCellID>* selection_)
                    : selection(selection_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                selection->push_back(i);
            }
        };

        Kokkos::parallel_for(this->getNbVec(), selectionPopulate(&conflictvertices));
        Kokkos::parallel_for(this->getNbVec(), ColorCell(this, &coloring, &conflictvertices));


        // get the conflicts
        struct ConflictSelection
        {
            const Graph *gr;
            Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic> > *coloring;
            kmds::GrowingView<kmds::TCellID> *selection;

            ConflictSelection(const Graph *gr_, Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic> > *coloring_, kmds::GrowingView<kmds::TCellID>* selection_)
                    : gr(gr_)
                    , coloring(coloring_)
                    , selection(selection_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
//                TCellID id = selection->get(i);
                TCellID id = i;

                Kokkos::View<kmds::TCellID *> neighbors;
                int nb = gr->getNeighbors(id, neighbors);

                int ownColor = (*coloring)[id];

                for(int ineighbor=0; ineighbor<nb; ineighbor++) {

                    TCellID neighborID = neighbors[ineighbor];

                    if((*coloring)[neighborID] == ownColor) {

                        // we keep only one of the two, the one with the highest id
                        if (neighborID < id) {
                            selection->push_back(id);
                            (*coloring)[id] = 0;
                            break;
                        }
                    }
                }
            }
        };

        // detect and solve the conflicts
        conflictvertices.clear();

        Kokkos::parallel_for(this->getNbVec(), ConflictSelection(this, &coloring, &conflictvertices));

        int nbConflicts = conflictvertices.getNbElems();

        while (nbConflicts > 0) {
            std::cout<<"nbConflicts "<<nbConflicts<<std::endl;

            Kokkos::parallel_for(nbConflicts, ColorCell(this, &coloring, &conflictvertices));

            conflictvertices.clear();

            Kokkos::parallel_for(this->getNbVec(), ConflictSelection(this, &coloring, &conflictvertices));
            nbConflicts =  conflictvertices.getNbElems();

            //exit(0);
        }

//        std::cout<<"graph::coloring "<<std::endl;
//        for(int i=0; i<this->getNbVec(); i++) {
//            std::cout<<i<<" "<<coloring[i]<<std::endl;
//        }

        struct ColorSelection
        {
            const Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic> > *coloring;
            std::vector<kmds::GrowingView<kmds::TCellID> >* color_sets;

            ColorSelection(const Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic> > *coloring_, std::vector<kmds::GrowingView<kmds::TCellID> >* color_sets_)
                    : coloring(coloring_)
                    , color_sets(color_sets_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                ((*color_sets)[(*coloring)[i]-1]).push_back(i);

            }
        };


        Kokkos::parallel_for(this->getNbVec(), ColorSelection(&coloring, &m_colorsets));

//        for(int i=0; i<nbMaxColors; i++) {
//            std::cout << "(*ASelections)[i] " << m_colorsets[i].getNbElems() << std::endl;
//
//            for(int j=0; j<m_colorsets[i].getNbElems(); j++) {
//                std::cout<<m_colorsets[i].get(j)<<" ";
//            }
//            std::cout<<std::endl;
//        }
    }
/*----------------------------------------------------------------------------*/
    kmds::GrowingView<kmds::TCellID>*
    Graph::getColoring(const int AColor)
    {
        return &m_colorsets[AColor];
    }
/*----------------------------------------------------------------------------*/
    void
    Graph::populateNonOriented(TFloat AMeanNbNeighbors) const
    {
        // random generator initialization
        //std::random_device rd;
        //std::mt19937 gen(rd);
        std::mt19937 gen(10);
        std::uniform_real_distribution<> dist(0., 1.);

        // neighbor threshold;
        double threshold = (double) AMeanNbNeighbors / (double) this->getNbVec();

        for(int id=0; id<this->getNbVec(); id++) {

            for(int ineighbor=id+1; ineighbor<this->getNbVec(); ineighbor++) {

                // vertex id has reached the maximum neighbor capacity
                if(m_nbneighbors(id) == m_maxneighbors) {
                    break;
                }
                // this potential neighbor vertex has reached the maximum neighbor capacity
                if(m_nbneighbors(ineighbor) == m_maxneighbors) {
                    continue;
                }

                double randomnumber = dist(gen);

                if(randomnumber < threshold) {

                    m_neighbors(id, m_nbneighbors(id)) = ineighbor;
                    m_nbneighbors(id)++;

                    m_neighbors(ineighbor, m_nbneighbors(ineighbor)) = id;
                    m_nbneighbors(ineighbor)++;
                }

            }

        }

        std::cout<<"end populate"<<std::endl;

    }
/*----------------------------------------------------------------------------*/
    void
    Graph::print() const
    {
        for(int id=0; id<this->getNbVec(); id++) {

            std::cout<<"id "<<id<<" | ";

            Kokkos::View<kmds::TCellID *> neighbors;
            int nb = this->getNeighbors(id, neighbors);

            for(int ineighbor=0; ineighbor<nb; ineighbor++) {
                std::cout<<neighbors[ineighbor]<<" ";
            }
            std::cout<<std::endl;
        }
    }
/*----------------------------------------------------------------------------*/
}  // namespace kmds
/*----------------------------------------------------------------------------*/
