/*----------------------------------------------------------------------------*/
#ifndef GMDS_PILLOW3D_H
#define GMDS_PILLOW3D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSSheet_export.h"
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  Pillow3D
     *  @brief  Class that provides a way to insert a full layer of hexes around
     *          a set of hexes, called a shrink set.
     *
     *          The mesh that is transformed must have R and N and the R2N
     *          fields. Geometric classification is not handled yet.
     */
    class GMDSSheet_API Pillow3D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         *  @param AMesh the mesh the pillowing operation will be performed on
         *  @param AEscapeBnd boolean value indicating the pillowing behaviour on the
         *                    boundary. True means that the pillow layer goes through
         *                    it. False means that the pillow layer goes along
         *                    boundary faces.
         */
        Pillow3D(Mesh* AMesh, const bool AEscapeBnd=true);
        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~Pillow3D();
        /*------------------------------------------------------------------------*/
        /** @brief Check if the mesh fits algorithm requirements, which are:
         *         - to be a 3D mesh,
         *         - to be a full hex mesh
         */
        bool isValid() const;
        /*------------------------------------------------------------------------*/
        /** @brief  Performs the pillowing
         * @param   ARegionIDs the set of regions we want to pillow
         * @param   AWithCheck if true check that the set of surfaces surrounding
         *          @ARegionIDs is valid
         *
         * @return true if the execution succeeded, false otherwise
         */
        bool execute(const std::vector<TCellID>& ARegionIDs, bool AWithCheck=false);
        /*------------------------------------------------------------------------*/
        /** @brief Performs the pillowing
         * @param  AFaces a set of virtual faces defined by nodes.
         * @param AWithCheck if true check that the set of surfaces @AFaces is valid
         *
         * @return true if the execution succeeded, false otherwise
         */
        bool execute(std::vector<VirtualFace>& AFaces, bool AWithCheck=false);
        /*------------------------------------------------------------------------*/
        /**@brief Setter of the boundary behaviour
         *
         * @param AEscape true means that the inserted sheet must go through the
         *                surface, while false means the sheet remains along the
         *                boundary surface.
         */
        void setBoundaryBehaviour(const bool AEscape);
        /*------------------------------------------------------------------------*/
        /**@brief Check if the surface enclose a shrink set
         * @param AFaces a set of faces defining a surface in the mesh
         * @return true if it encloses, false otherwise
         */
        bool checkSurface(std::vector<VirtualFace>& AFaces);
    private:
        /*------------------------------------------------------------------------*/
        /**@brief   Orient faces of @AFaces in a consistent manner
         *
         *          Does not require the N2R connectivity, only works on the faces
         *
         * @param AFaces
         */
        void orient(std::vector<VirtualFace>& AFaces);
        /*------------------------------------------------------------------------*/
        /** @brief Fill the local connectivity N2R
         */
        void buildLocalN2R();
        /*------------------------------------------------------------------------*/
        /** @brief Check if a node is on the mesh boundary
         * @param ANodeID the id of the node to be checked
         */
         bool isABoundaryNode(const TCellID ANodeID);
        /*------------------------------------------------------------------------*/
        /** @brief Check if a quad face is on the mesh boundary
         * @param AN1 id of the first face node
         * @param AN2 id of the second face node
         * @param AN3 id of the third face node
         * @param AN4 id of the fourth face node
         */
        bool isABoundaryFace(const TCellID AN1,const TCellID AN2,
                             const TCellID AN3,const TCellID AN4);
        /*------------------------------------------------------------------------*/
        /** @brief Check if a quad face is on the mesh boundary
         * @param AF A virtual  face
         */
        bool isABoundaryFace(const VirtualFace& AF);

        /*------------------------------------------------------------------------*/
        /** @brief Get the regions sharing face AF
         * @param AF a virtual face
         * @return the adjacent regions (1 or 2)
         */
        std::vector<TCellID> getAdjacentRegions(const VirtualFace& AF);

        /*------------------------------------------------------------------------*/
        /** @brief Get the regions sharing a face with region AID
         * @param AID a region id
         * @return the adjacent regions (1 or 2)
         */
        std::vector<TCellID> getRegionsSharingAFace(const TCellID ARegionID);


        /*------------------------------------------------------------------------*/
        /** @brief Get the regions sharing nodes @AN1, @AN2, @AN3 and @AN4
         * @param AN1 id of the first face node
         * @param AN2 id of the second face node
         * @param AN3 id of the third face node
         * @param AN4 id of the fourth face node
         */
        std::vector<TCellID> getAdjacentRegions(const TCellID AN1,const TCellID AN2,
                                                const TCellID AN3,const TCellID AN4);
    private:
        /** a mesh */
        Mesh* m_mesh;
        /** options to indicate the behaviour of the operation along the boundary*/
        bool m_escape_bnd;
        /** Inverse connectivity N23 build on purpose for this algorithm*/
        std::map<TCellID, std::vector<TCellID> > m_N2R;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_PILLOW3D_H
