/*----------------------------------------------------------------------------*/
/*
 * Connectivity.h
 *
 *  Created on: 03/13/2017
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_CONNECTIVITY_H_
#define KMDS_CONNECTIVITY_H_
/*----------------------------------------------------------------------------*/
// STL headers
#include <string>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
class Mesh;
/*----------------------------------------------------------------------------*/
class EXPORT_KMDS Connectivity
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor
         *
         * \param[in] AOwner the mesh that owns this cell container
         * \param[in] AType  the type of connectivity
         * \param[in] ACapacity initial capacity of the container
         */
        Connectivity(Mesh* AOwner, const EMeshDefinition AT, const TCellID ACapacity = 16);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~Connectivity();

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the connectivity definition
         */
        EMeshDefinition getType() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the connectivity name
         */
        std::string getName() const;

        /*------------------------------------------------------------------------*/
        /** \brief Returns the connectivity capacity in terms of single TCellID
         *         items
         */
        TSize getAbsoluteCapacity() const;

        /*------------------------------------------------------------------------*/
        /** \brief Returns the connectivity capacity in terms of number of records
         */
        TCellID getCapacity() const;

        /*------------------------------------------------------------------------*/
        /** \brief Returns the connectivity capacity in terms of single TCellID
         *         items
         */
        void setAbsoluteCapacity(const TSize ASize);

        /*------------------------------------------------------------------------*/
        /** \brief Returns the connectivity capacity in terms of number of records;
         *         must be higher than the highest id.
         */
        void setCapacity(const TCellID ASize);

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the connectivity name associated to a type
         */
        static std::string getName(const EMeshDefinition AD);

        void set(const TCellID AKey, const Kokkos::View<TCellID*>& AV);
        void set(const TCellID AKey, const TCellID* AV, const TSize ASize);

        void get(const TCellID AKey, Kokkos::View<TCellID*>& AV) const;

        void reinit();

        void setTop(TSize ASize);
        void setTop_C2C(TSize ASize);
        Kokkos::View<TSize*>& getFirst();
        Kokkos::View<int*>& getNbCells();
        Kokkos::View<TCellID*>& getC2C();


 private:
        /** mesh owner of this cell container */
        Mesh* m_mesh;
        /** type of connectivity*/
        EMeshDefinition m_type;

        /** First connected cell, if = NullID means empty*/
        Kokkos::View<TSize*> m_first;
        /** Nb connected cell */
        Kokkos::View<int*> m_nb_cells;
        /** Stored connectivities */
        Kokkos::View<TCellID*> m_C2C;
        /** top of the heap including holes so */
        TSize m_top;
        TSize m_top_C2C;
};
/*----------------------------------------------------------------------------*/
class EXPORT_KMDS ConnectivityHelper
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor
         *
         * \param[in] AOwner the mesh that owns this cell container
         */
        ConnectivityHelper(Mesh* AOwner);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~ConnectivityHelper();

        /*------------------------------------------------------------------------*/
        /** \brief Build the connectivity N2E for all the nodes and edges available
         *         in m_mesh. Must not be done in a parallel kernel
         */
        void buildN2E();

        /*------------------------------------------------------------------------*/
        /** \brief Build the connectivity N2F for all the nodes and faces available
         *         in m_mesh. Must not be done in a parallel kernel
         */
        void buildN2F();

        /*------------------------------------------------------------------------*/
        /** \brief Build the connectivity N2R for all the nodes and regions available
         *         in m_mesh. Must not be done in a parallel kernel
         */
        void buildN2R();
        void buildN2R_variant_0();
        void buildN2R_variant_1();

        /*------------------------------------------------------------------------*/
        /** \brief Build the edges and connectivity E2F for all the faces available
         *         in m_mesh. Must not be done in a parallel kernel
         */
        // TODO sequential for now
        void buildEandE2F();
        void buildEandE2F_2D_variant_0();
        void buildEandE2F_2D_variant_1();

        // this variant does not assume a max of two faces per edges
        void buildEandE2F_3D();

    /*------------------------------------------------------------------------*/
    /** \brief Build the edges and connectivity E2R for all the regions available
     *         in m_mesh. Must not be done in a parallel kernel
     */
    void buildEandE2R();

        /*------------------------------------------------------------------------*/
        /** \brief Build the faces and connectivity F2R for all the regions available
         *         in m_mesh. Must not be done in a parallel kernel
         */
         // TODO sequential for now
        void buildFandF2R();
        void buildFandF2R_variant_0();

        /*------------------------------------------------------------------------*/
        /** \brief Build the connectivity F2F_byN for all the faces available
         *         in m_mesh. Must not be done in a parallel kernel
         */
        void buildF2F_byN(const kmds::Connectivity* c_N2F);

        /*------------------------------------------------------------------------*/
        /** \brief Build the connectivity F2F (adjacency by E) for all the faces available
        *         in m_mesh. Must not be done in a parallel kernel
        */
        void buildF2F(const kmds::Connectivity* c_N2F);

        /*------------------------------------------------------------------------*/
        /** \brief Build the connectivity R2R_byN for all the regions available
         *         in m_mesh. Must not be done in a parallel kernel
         */
        void buildR2R_byN(const kmds::Connectivity* c_N2R);

        /*------------------------------------------------------------------------*/
        /** \brief Build the connectivity N2N for all the nodes available
         *         in m_mesh. Based on N2E2N adjacency, with edges found in ACellType
         *         Must not be done in a parallel kernel
         */
        // TODO sequential for now
        void buildN2N(const EMeshDefinition ACellType);


 private:
        /** mesh owner of the connectivity */
        Mesh* m_mesh;
};
/*----------------------------------------------------------------------------*/
}  // namespace KMDS
/*----------------------------------------------------------------------------*/
#endif /* KMDS_CONNECTIVITY_H_ */
/*----------------------------------------------------------------------------*/
