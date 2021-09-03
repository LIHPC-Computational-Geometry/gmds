/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
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
/*----------------------------------------------------------------------------*/
/** \file    FakeTypes.h
 *  \author  legoff
 *  \date    02/23/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_FAKETYPES_H_
#define KMDS_FAKETYPES_H_
/*----------------------------------------------------------------------------*/
// STL headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Kokkos headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// KMDS headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
    /*----------------------------------------------------------------------------*/
    /** \class FakeEdge
     *
     *  \brief This class stores an edge relatively to its node ids.
     *  	   The node order is given starting from the lowest node id, N1.
     *
     */
    class FakeEdge{

            public:
            /*-----------------------------------------------------------------------*/
            /** \class ID
             *
             *  \brief defines an id for a fake edge
             */
            struct EdgeID{

                EdgeID(const TCellID& AID1=1, const TCellID& AID2=1)
                { id_tab_[0]=AID1; id_tab_[1]=AID2;}

                void set(const TCellID& AID1=1, const TCellID& AID2=1)
                { id_tab_[0]=AID1; id_tab_[1]=AID2;}

                bool operator==(const EdgeID& AID) const{
                    return (id_tab_[0]==AID.id_tab_[0] && id_tab_[1]==AID.id_tab_[1] );
                }

                bool operator<(const EdgeID& AID) const{
                    return ((id_tab_[0]<AID.id_tab_[0]) ||
                            (id_tab_[0]==AID.id_tab_[0] &&
                                    id_tab_[1]<AID.id_tab_[1]) );
                }

                TCellID getID1() const {return id_tab_[0];}
                TCellID getID2() const {return id_tab_[1];}

                private:
                TCellID id_tab_[2];
            };

            FakeEdge(){
                init(NullID,NullID);
            }

            FakeEdge(const TCellID& AID1, const TCellID& AID2){
                init(AID1,AID2);
            }

            ~FakeEdge(){;}

            std::vector<TCellID> node_ids() const {
                std::vector<TCellID> node_ids;
                node_ids.push_back(id_.getID1());
                node_ids.push_back(id_.getID2());
                return node_ids;
            }

            TCellID first() const {return id_.getID1();}
            TCellID second() const {return id_.getID2();}
            void set(const TCellID& AID1, const TCellID& AID2){
                init(AID1,AID2);
            }

            EdgeID getID() const {return id_;}

            bool operator==(const FakeEdge& AF)
            {	return (id_==AF.id_);}

            bool operator!=(const FakeEdge& AF)
            {	return !(id_==AF.id_);}

            bool operator<(const FakeEdge& AF) const
            {	return (id_<AF.id_);}

            private:

            void init(const TCellID& AID1, const TCellID& AID2){
                if(AID1<AID2)
                    id_.set(AID1,AID2);
                else
                    id_.set(AID2,AID1);

            }


            private:

            EdgeID id_;
    };
    /*----------------------------------------------------------------------------*/
    /** \class FakeFace
     *
     *  \brief This class stores the node ids of a face. The node order is given
     *		   starting from the lowest node id, N1,  and traversing the nodes
     *		   in the direction of the N1 incident node having the lowest id.
     *
     */
    class FakeFace{

            public:


            /*-----------------------------------------------------------------------*/
            /** \class ID
             *
             *  \brief defines an id made of 3 id
             */
            struct FaceID{

                FaceID(){;}
                FaceID(const std::vector<TCellID>& AIDs){tab_ = AIDs;}

                void set(const std::vector<TCellID>& AIDs){tab_ = AIDs;}

                bool operator==(const FaceID& AID) const{
                    if (tab_.size()!=AID.tab_.size())
                        return false;

                    for(unsigned int i=0;i<tab_.size();i++)
                        if(tab_[i]!=AID.tab_[i])
                            return false;

                    return true;
                }

                bool operator<(const FaceID& AID) const{
                    if(tab_.size()<AID.tab_.size())	return true;
                    if(tab_.size()>AID.tab_.size()) return false;
                    /* now same size */
                    for(unsigned int i=0;i<tab_.size();i++)
                        if(tab_[i]<AID.tab_[i])		  return true;
                        else if (tab_[i]>AID.tab_[i]) return false;
                    /* this and AID are then equals */
                    return false;
                }

                std::vector<TCellID> getNodeIDs() const {return tab_;}

                private:
                std::vector<TCellID> tab_;
            };

            FakeFace(const std::vector<TCellID>& ANodeIDs){
                init(ANodeIDs);
            }

            ~FakeFace(){;}

            FaceID getID() const {return id_;}

            std::vector<TCellID> node_ids() const {
                return id_.getNodeIDs();
            }


            void set(const std::vector<TCellID>& ANodeIDs){
                init(ANodeIDs);
            }

            bool operator==(const FakeFace& AF) const
            {	return (id_==AF.id_);}

            bool operator<(const FakeFace& AF) const
            {	return (id_<AF.id_);}

            std::vector<FakeEdge> getFakeEdges() const {
                std::vector<FakeEdge> fakeEdges;
                std::vector<TCellID> nids = this->node_ids();
                for(int i=0; i<nids.size(); i++) {
                    fakeEdges.push_back(FakeEdge(nids[i],nids[(i+1)%nids.size()]));
                }
                return fakeEdges;
            }

            std::vector<TCellID> getTriangle(TCellID AID) const {
                std::vector<TCellID> tri(3);
                tri[0] = AID;

                std::vector<TCellID> ids = this->node_ids();
                for(int i=0; i<ids.size(); i++) {
                    if(AID == ids[i]) {
                        tri[1] = ids[(i+1)%ids.size()];
                        tri[2] = ids[(i-1 + ids.size())%ids.size()];
                    }
                }
                return tri;
            }

            bool hasNode(TCellID AID) const {
                std::vector<TCellID> nids = this->node_ids();
                for(auto n: nids) {
                    if(AID == n) {
                        return true;
                    }
                }
                return false;
            }


            private:

            void init_front(const std::vector<TCellID>& ANodeIDs, const TSize & AFirstItem,
            std::vector<TCellID>& AOutNodeIDs){
                for(unsigned int i = AFirstItem;i<ANodeIDs.size();i++)
                    AOutNodeIDs.push_back(ANodeIDs[i]);
                for(int i = 0;i<AFirstItem;i++)
                    AOutNodeIDs.push_back(ANodeIDs[i]);
            }
            void init_back(const std::vector<TCellID>& ANodeIDs, const TSize& AFirstItem,
            std::vector<TCellID>& AOutNodeIDs){
                for(int i = AFirstItem;i>=0;i--)
                    AOutNodeIDs.push_back(ANodeIDs[i]);
                for(int i = ANodeIDs.size()-1;i>AFirstItem;i--)
                    AOutNodeIDs.push_back(ANodeIDs[i]);
            }
            void init(const std::vector<TCellID>& ANodeIDs){
                /* ANodeIDs gives an ordered description of a face by nodes */
                std::vector<TCellID> node_ids;
                node_ids.reserve(ANodeIDs.size());

                /* look for the smallest node id*/
                TCellID lowest_id= NullID;
                TSize     lowest_item = 0;
                TSize nb_ids = ANodeIDs.size();
                for(int i=0;i<nb_ids;i++)
                    if (ANodeIDs[i]<lowest_id)
                    {
                        lowest_id = ANodeIDs[i];
                        lowest_item = i;
                    }
                /* study of the surrounding node*/
                TCellID left_value, right_value;
                if(lowest_item==0){
                    left_value = ANodeIDs[nb_ids-1];
                    right_value = ANodeIDs[1];
                }
                else if(lowest_item==nb_ids-1){
                    left_value = ANodeIDs[nb_ids-2];
                    right_value = ANodeIDs[0];
                }
                else{
                    left_value = ANodeIDs[lowest_item-1];
                    right_value = ANodeIDs[lowest_item+1];
                }
                if(left_value<right_value)
                    init_back(ANodeIDs,lowest_item,node_ids);
                else
                    init_front(ANodeIDs,lowest_item,node_ids);

                id_.set(node_ids);

            }


            private:

            FaceID id_;
    };
    /*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* KMDS_FAKETYPES_H_ */
/*----------------------------------------------------------------------------*/
