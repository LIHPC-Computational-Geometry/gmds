/*----------------------------------------------------------------------------*/
#ifndef POLYCUBE_TOOLBOX_H
#define POLYCUBE_TOOLBOX_H
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/ig/Mesh.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
/*----------------------------------------------------------------------------*/
using namespace std;
/*----------------------------------------------------------------------------*/
enum labelling {Xp,Yp,Zp,Xm,Ym,Zm};
/*----------------------------------------------------------------------------*/
Eigen::Vector3d math2eigen(gmds::math::Vector v);

const gmds::math::Vector vect[6] = {gmds::math::Vector(1,0,0),
                                    gmds::math::Vector(-1,0,0),
                                    gmds::math::Vector(0,1,0),
                                    gmds::math::Vector(0,-1,0),
                                    gmds::math::Vector(0,0,1),
                                    gmds::math::Vector(0,0,-1) };
/*----------------------------------------------------------------------------*/
namespace gmds {

    class PolycubeToolbox {
    private:
        Mesh *mesh;
        int mark_smooth_bnd_nodes;
        int mark_bnd_nodes;
        int mark_bnd_faces;
        int mark_sharp_edges;
        int mark_bnd_edges;

        double BOUNDING_BOX_LENTH;

        gmds::Variable<math::Vector> *face_normals;
        gmds::Variable<int> *face_chosenNormal_ID;
        gmds::Variable<math::Vector> *face_chosenNormal;
        gmds::Variable<int> *face_chart;
        gmds::Variable<int> *face_constrained;

        gmds::Variable<math::Vector> *point_normals;
        gmds::Variable<math::Vector> *point_chosenNormal;
        gmds::Variable<Eigen::Matrix<double, 3, 3>> *point_rotationMatrix;
        gmds::Variable<Eigen::Quaternion<double>> *point_quaternion;

        gmds::Variable<int> *node_updatedX;
        gmds::Variable<int> *node_updatedY;
        gmds::Variable<int> *node_updatedZ;

        gmds::Variable<int> *node_chartCorner;
        gmds::Variable<int> *node_turningPoint;
        gmds::Variable<int> *node_browsed;

        //this variable is of use to compute turning point
        gmds::Variable<math::Vector> *edge_normal;


        double ANGLE_DETERMINENT;

        typedef std::vector<TCellID> chart;
        typedef std::vector<TCellID> chart_edge;
        struct Triplet_charts_edge {
            chart_edge edge;
            int chart1, chart2;
            TCellID corner1, corner2;
            std::vector<TCellID> turningPoints;
        };

        std::vector<chart> charts;
        std::vector<Triplet_charts_edge> charts_edges;


        //--------------Gregson2011------------------//
        Eigen::SparseLU<Eigen::SparseMatrix<double>> poisson_solver;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> adjacency_solver;
        Eigen::SparseMatrix<double> poisson_matrix;
        Eigen::SparseMatrix<double> adjacency_matrix;
        int GREGSON2011_ITERATION_NUMBER;
        bool GREGSON2011_WRITING_RESULT;
        bool GREGSON2011_IS_INIT;
        bool GREGSON2011_RECOMPUTING_VARIABLES_IN_ITERATIONS;
        //------------------------------------------//


        std::vector<TCellID> bnd_faces;
        map<TCellID, int> linker_ID2bnd_faces; // useful to know to which id was assigned a face.
        std::vector<bool> gc_face_is_constrained;
        std::map<std::pair<TCellID, int>, double> gc_face_constraint;


    public:
        PolycubeToolbox(Mesh *mesh);

        void transform_mesh_to_polycube();

        Mesh getMesh();


        void set_angle_determinent(double);

        // TODO rotate the mesh instead of this
        // void set_1st_axis(math::Vector);

        void gregson2011Init();

        void gregson2011SetWritingResult(bool);

        void gregson2011SetVariableComputationInIteration(bool);

        void gregson2011Iteration();

        void gregson2011VariablesGeneration();

        void gregson2011Run(int nb_iteration = 5);


        void smoothing();


        void chartsGeneration(bool speaking = true);

        void compute_turning_point(bool speaking = true);

        void projectionPreprocess();

        void projectOnPolycube();

        bool checkingChartCornerOverloading();

        bool checkingChartDivision();

        void graphCutBasedFaceSegmentation();

        void
        graphCutHillClimbing(double ADeltaStart = 0, double ADeltaIncr = 0.01, double bounding_box_ratio = 0.05);

        void graphCutRun();

        void graph_cut_labelling_update();

        void computeNodeInformation();

        void test_fonction();

        /**
            This function return a vector containing the point's TCellIDs of the
            corner of each faces of the blocs. The TCellID don't change during the
            algorithm. The bloc on a reloaded mesh will still be valid.
        */
        std::vector<std::vector<TCellID>> get_faces_of_block();

    private:
        void markOnBoundary();

        void updateBoundary();

        Eigen::SparseMatrix<double> computeAdjacency();

        Eigen::SparseMatrix<double> computePoisson();

        Eigen::SparseMatrix<double> computeLaplace();

        labelling closestNormal(math::Vector v);


        static double graphCutCpq(int p, int q, int label_p, int label_q, void *extraData);

        inline math::Vector convert_edge_to_vector(Edge);

        inline double angleFace2Node(TCellID face, TCellID node);

    };
}
/*----------------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------------*/
