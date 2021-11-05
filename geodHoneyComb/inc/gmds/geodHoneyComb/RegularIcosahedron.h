/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 04/11/2021.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_REGULARICOSAHEDRON_H
#define GMDS_REGULARICOSAHEDRON_H
/*----------------------------------------------------------------------------*/
#include <memory>
#include <gmds/ig/Mesh.h>
#include "LIB_GMDS_GEOD_HONEY_COMB_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    /** @class A regular icosahedron is a convex polyhedron with 20 faces, 30
     *          edges and 12 vertices. It is one of the five Platonic solids,
     *          and the one with the most faces.
     *          We build it as a surface mesh made of nodes and triangles and the
     *          N2F and F2N connections.
     */
    class LIB_GMDS_GEOD_HONEY_COMB_API RegularIcosahedron{
    public:
        /** @brief Build a regular icosahedron with radius \p ARadius and \p ACenter
         *
         * @param ACenter polyhedron center
         * @param ARadius radius of the polyhedron circumsphere.
         * @param ANBSubdivision the number of time we subdivide the icosahedron
         */
        RegularIcosahedron(const math::Point& ACenter=math::Point(0,0,0),
                           const TCoord ARadius=1,
                           const int ANBSubdivision=1);

        /** @brief Destructor
         */
        virtual ~RegularIcosahedron();

        /** @brief Create the dual icosahedron with each face splitted in quads
         */
        void performQuadDualization();

        /** @brief Gives access to the mesh representation of the icosahedron
         */
        std::unique_ptr<Mesh> getRepresentation();

    private:
        /** @brief position a point on the circumsphere of the polyhedron
         *  @param  APoint point to project
         *  @return resulting point
         */
        math::Point projectOnSphere(const math::Point& APoint);
        /** @subdivide the original icosahedron according to m_subdivision value
         */
        void subdivide();
    private:
        /** the number of times we subdivide the original icosahedron */
        int m_subdivision;
        /** center of the circumsphere*/
        math::Point m_center;
        /** radius of the circumsphere*/
        TCoord m_radius;
        /** Surface of the polyhedron represented as a mesh*/
        std::unique_ptr<Mesh> m_representation;

    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_REGULARICOSAHEDRON_H
/*----------------------------------------------------------------------------*/
