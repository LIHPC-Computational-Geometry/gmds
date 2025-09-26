/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEODHEXMESHER_H
#define GMDS_GEODHEXMESHER_H
/*----------------------------------------------------------------------------*/
#include <map>
#include <memory>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/DiscretizationScheme1D.h>
#include "GMDSGeodHoneyComb_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    /** @class This class gathers some basic algorithm to build a honeycomb
     *          mesh of a spherical domain
     */
    class GMDSGeodHoneyComb_API GeodHexMesher{
    public:
        typedef enum {
            GEOD_SUCCESS,
            GEOD_FAILURE,
            GEOD_FAILURE_WRONG_RADIUS,
            GEOD_FAILURE_LAYER_RADIUS_NEGATIVE,
            GEOD_FAILURE_LAYER_RADIUS_GREATER_THAN_SPHERE_RADIUS
        } OpResult;

        typedef enum {
            GEOD_NOT_DEFINED,
            GEOD_SHELL,
            GEOD_SPHERE
        } geodType;


        /** @brief Constructor
         **/
        GeodHexMesher();

        /** @brief Destructor
         */
        virtual ~GeodHexMesher();

        void setRadius(const double AR);
        void setCenter(const math::Point& AP);

        bool checkValidity() const;

        /** @brief set the layer data info.
         *
         * Layers radius must be enclosed between 0 and outer radius, both excluded
         * If provided layers are not give as expected, an exception is raised.
         * If the smallest radius is 0, we will generate a sphere, otherwise a shell.
         * the discretization associated to a radius corresponds to the discretization scheme
         * applied between this radius value to the next upper radius value.
         *
         * The highest layer radius must be < m_radius
         *
         * @param ALayerData couples of radius, discretization data
         */
         void setLayerData(std::map<double, math::DiscretizationScheme1D*>& ALayerData);
        /*-------------------------------------------------------------------*/
        /** @brief Execute the algorithm to get the final mesh
         * @return The result of the operation (success or failure)
         */
        OpResult execute();


        /** @brief Gives access to the final mesh representation. */
        std::unique_ptr<Mesh> getMesh();


    private:

        /** sphere radius, or exter radius of a shell */
        double m_radius;
        /** sphere center*/
        math::Point m_center;
        /** we generate a sphere or a shell depending on this value*/
        bool m_is_sphere;

        /** strict ordered list of radius values from the center to the
         *  external of the sphere. The minimal values must be > m_inner_radius
         *  and the maximal one must be < to m_radius. For each zone a discretization
         *  data is required too.
         */
        std::map<double,math::DiscretizationScheme1D*> m_layer_data;

        /** obtained mesh*/
        std::unique_ptr<Mesh> m_mesh;
    };
    /*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_GEODHEXMESHER_H
/*----------------------------------------------------------------------------*/
