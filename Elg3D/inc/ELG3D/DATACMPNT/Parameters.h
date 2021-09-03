/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Parameters.h
 *  \author  legoff
 *  \date    04/18/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_PARAMETERS_H_
#define ELG3D_PARAMETERS_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
// elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
    typedef enum {
        glpk,
        r3d,
        ortho
    } Parameter_badpillowing_method;

/*----------------------------------------------------------------------------*/
struct Parameters{

    // set the dimension
    static int dimension;

    // set the cell type (will be KMDS_FACE or KMDS_REGION)
    static kmds::ECellType celltype;

    // set if the bad pillowing configurations are to be avoided (valid only in 3D)
    static bool badpillowing_detection_ON;

    // set if the bad pillowing configurations are to be avoided (valid only in 3D)
    static Parameter_badpillowing_method badpillowing_detection_method;

    // the distance to the seed nodes cells will be marked for pillowing (currently only in 2D)
    static int pillow_node_seed_range;

    // the quality criteria below which some of the steps will not go
    static double quality_threshold;

    // prefix name to the outputs of the code
    static std::string output_prefix;


    // DEBUG PARAMETERS
    static int debug_verbosity_lvl;

    static bool debug_ouputfile_pillow;

    static bool debug_ouputfile_pixels;

    // mapping between pair of materials indices (lower, higher) to an index
    // used for debug, outputs purposes
    static std::map<std::pair<int, int>, int> mat_pair_2_int;
};

/*----------------------------------------------------------------------------*/
// 3-Refinement tables
    template<int TSize>
    class tableMarkedNodes {
    public :
        tableMarkedNodes() {
            for(unsigned int iMark = 0; iMark<TSize; iMark++) {
                mark_[iMark] = 0;
            }
        };

        tableMarkedNodes(const bool& AMark1, const bool& AMark2, const bool& AMark3, const bool& AMark4) {

            mark_[0] = AMark1;
            mark_[1] = AMark2;
            mark_[2] = AMark3;
            mark_[3] = AMark4;
        };

        tableMarkedNodes(const bool& AMark1, const bool& AMark2, const bool& AMark3, const bool& AMark4,
                         const bool& AMark5, const bool& AMark6, const bool& AMark7, const bool& AMark8) {

            mark_[0] = AMark1;
            mark_[1] = AMark2;
            mark_[2] = AMark3;
            mark_[3] = AMark4;
            mark_[4] = AMark5;
            mark_[5] = AMark6;
            mark_[6] = AMark7;
            mark_[7] = AMark8;
        };

        tableMarkedNodes(const bool AMarks[TSize]) {
            for(unsigned int iMark = 0; iMark<TSize; iMark++) {
                mark_[iMark] = AMarks[iMark];
            }
        };

        tableMarkedNodes(const tableMarkedNodes<TSize>& ATableMarkedNodes) {
            for(unsigned int iMark=0; iMark<TSize; iMark++) {
                mark_[iMark] = ATableMarkedNodes.mark_[iMark];
            }
        };

        ~tableMarkedNodes() {

        };

        tableMarkedNodes<TSize>& operator=(const tableMarkedNodes<TSize>& ATableMarkedNodes) {
            for(unsigned int iMark=0; iMark<TSize; iMark++) {
                mark_[iMark] = ATableMarkedNodes.mark_[iMark];
            }

            return *this;
        };

        bool operator<(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
            for(unsigned int iMark=0; iMark<TSize; iMark++) {
                if(!mark_[iMark] && ATableMarkedNodes.mark_[iMark]) {
                    return true;
                }
                if(mark_[iMark] && !ATableMarkedNodes.mark_[iMark]) {
                    return false;
                }
            }

            return false;
        };

        bool operator>(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
            for(unsigned int iMark=0; iMark<TSize; iMark++) {
                if(mark_[iMark] && !ATableMarkedNodes.mark_[iMark]) {
                    return true;
                }
                if(!mark_[iMark] && ATableMarkedNodes.mark_[iMark]) {
                    return false;
                }
            }

            return false;
        };

        bool operator!=(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
            for(unsigned int iMark=0; iMark<TSize; iMark++) {
                if(mark_[iMark] != ATableMarkedNodes.mark_[iMark]) {
                    return true;
                }
            }

            return false;
        };

        bool operator==(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
            for(unsigned int iMark=0; iMark<TSize; iMark++) {
                if(mark_[iMark] != ATableMarkedNodes.mark_[iMark]) {
                    return false;
                }
            }

            return true;
        };

        bool isMarked(const int& AIndex) const {
            return mark_[AIndex];
        };

        bool isEqual(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
            return (*this==ATableMarkedNodes);
        };

        bool get(const int AIndex) const {
            return mark_[AIndex];
        }

    private :
        bool mark_[TSize];
    };

    template<int TSize>
    class compare_tableMarkedNodes {
    public :
        bool operator()(const tableMarkedNodes<TSize>& ATableMarkedNodes1, const tableMarkedNodes<TSize>& ATableMarkedNodes2) {
            return (ATableMarkedNodes1<ATableMarkedNodes2);
        }
    };
    struct Refinement_tables {
        static std::set<tableMarkedNodes<4>, compare_tableMarkedNodes<4> > lookupNodesValid_3_2D;
        static std::map<tableMarkedNodes<4>, tableMarkedNodes<4>, compare_tableMarkedNodes<4> > lookupNodesSolve_3_2D;
        static std::set<tableMarkedNodes<8>, compare_tableMarkedNodes<8> > lookupNodesValid_3_3D;
        static std::map<tableMarkedNodes<8>, tableMarkedNodes<8>, compare_tableMarkedNodes<8> > lookupNodesSolve_3_3D;

        static std::array<int, 8> nodesOnFaceFromExternal_3ab_2D;
        static std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > nbNodesOnFaceToBuild_3a_2D;
        static std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > nbFacesOnFaceToBuild_3a_2D;
        static std::map<tableMarkedNodes<4>, std::vector<std::array<int, 2> >, compare_tableMarkedNodes<4> > nodesOnFaceToBuild_3a_2D;
        static std::map<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > facesOnFaceToBuild_3a_2D;
        static std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > nbNodesOnFaceToBuild_3b_2D;
        static std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > nbFacesOnFaceToBuild_3b_2D;
        static std::map<tableMarkedNodes<4>, std::vector<std::array<int, 2> >, compare_tableMarkedNodes<4> > nodesOnFaceToBuild_3b_2D;
        static std::map<tableMarkedNodes<4>, std::vector< std::array<int, 4> > > facesOnFaceToBuild_3b_2D;
        static std::map<tableMarkedNodes<4>, int, compare_tableMarkedNodes<4> > nbNodesOnFaceToBuild_3b_3D;
        static std::array<int, 48> nodesOnRegionFromExternal_3b_3D;
        static std::map<tableMarkedNodes<8>, int, compare_tableMarkedNodes<8> > nbNodesOnRegionFromExternal_3b_3D;
        static std::map<tableMarkedNodes<8>, int, compare_tableMarkedNodes<8> > nbNodesOnRegionToBuild_3b_3D;
        static std::map<tableMarkedNodes<8>, std::vector<std::array<int, 3> >, compare_tableMarkedNodes<8> > nodesOnRegionToBuild_3b_3D;
        static std::map<tableMarkedNodes<8>, int, compare_tableMarkedNodes<8> > nbRegionsOnRegionToBuild_3b_3D;
        static std::map<tableMarkedNodes<8>, std::vector< std::array<int, 8> > > regionsOnRegionToBuild_3b_3D;
        static std::map<tableMarkedNodes<4>, std::array<int, 12> > nodesOnRegionToBuild_3a_3D;
        static std::map<tableMarkedNodes<8>, std::array<int, 8>, compare_tableMarkedNodes<8> > permutNodesOfRegion_3b_3D;
    };
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_PARAMETERS_H_ */
