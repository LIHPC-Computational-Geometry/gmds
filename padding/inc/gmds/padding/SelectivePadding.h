/*----------------------------------------------------------------------------*/
#ifndef SELECTIVE_PADDING_H_
#define SELECTIVE_PADDING_H_
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace  gmds{
/*----------------------------------------------------------------------------*/
/** \class  SelectivePadding
 *  \brief  ???
 */
    class  SelectivePadding
    {
    public:

        enum Option{
            SBP        //SIMPLE BINARY PROBLEM
        };
        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *  @param AMesh the mesh where we work on
         */
        SelectivePadding(gmds::Mesh* AMesh);

        /*------------------------------------------------------------------------*/
        /** @brief Set the variable given hard faces as an input.
         *  @param AVar a binary face variable: 1 means the face is hard
         *              constrained, 0 means not constrained.
         */
        void setHardFaces(Variable<int>* AVar);
        /*------------------------------------------------------------------------*/
        /** @brief Set the variable given the faces to be padded as an output.
         *  @param AVar a binary face variable: 1 means the face has to be padded,
         *              0 means no.
         */
        void setPaddingFaces(Variable<int>* AVar);

        /*------------------------------------------------------------------------*/
        /** @brief Function to be called for running the face selectio process
         *  @param AOption options to define the way of selection (SBP bu default)
         */
        void execute(const Option AOption=Option::SBP);

        /**@brief Set the flag indicating to write (true) or not (false) the
         *        problem to be solved
         * @param AWithOutput Boolean flag
         */
        void setCplexOutputFlag(const bool AWithOutput=false);

        /**@brief Set the file name where writing the problem to be solved
         * @param AFileName file name
         */
        void setCplexFileName(const std::string& AFileName="prb.lp");
    private:

        /*------------------------------------------------------------------------*/
        /** \brief Contains the problem construction and execution using GLPK.
         */
        void buildAndSolveSBP();
    private:

        /** Background triangular mesh */
        gmds::Mesh* m_mesh;
        /** hard constraint on faces */
        Variable<int>* m_hard_constraint;
        /** faces to be padded */
        Variable<int>* m_padding;

        /** flag to write the solved problem in cplex format */
        bool m_cplex_with_output;
        /** name of the cplex output file*/
        std::string m_cplex_file_name;

    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //SELECTIVE_PADDING_H_
/*----------------------------------------------------------------------------*/
