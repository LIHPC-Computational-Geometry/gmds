/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/8/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ABSTRACT_SERVICE_H
#define GMDS_ABSTRACT_SERVICE_H
/*----------------------------------------------------------------------------*/
#include <set>
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds {
    class Property;
    class AbstractData;
    /*------------------------------------------------------------------------*/
    /** @class AbstractService
     *
     * @brief Defines the minimal interface of each service implemented in
     *        GMDS
     *
     */
    class AbstractService {
    public:

        /*--------------------------------------------------------------------*/
        /** @brief Launch the service execution without doing any check
         */
        virtual void execute()=0;


        /*--------------------------------------------------------------------*/
        /** @brief Launch the service execution after checking input data
         *
         * @throw a GMDSException if the check fails
         */
        void executeAfterChecking();

        /*--------------------------------------------------------------------*/
        /** @brief Add a data input
         * @param AData
         */
        void addInput(const AbstractData* AData);

        /**@brief Add a constraint on @AData, which is that @AData must
         *        satisfy the property @AProp
         *
         * @param AData data we want to put a constraint on
         * @param AProp the constraint definition
         *
         * @return true if the constraint is added, false if @AData is not
         *              an input or output of the service
         */
        bool addConstraint(const AbstractData* AData,
                           const Property* AProp);

        /**@brief check the validity of input data
         *
         * @return true if input data are okay, false otherwise
         */
        bool checkInput();

        /**@brief check the validity of output data
         *
         * @return true if output data are okay, false otherwise
         */
        bool checkOutput();



    private:

        /** input data of the service*/
        std::set<const AbstractData*> m_input;

        /** output data of the service*/
        std::set<const AbstractData*> m_output;

        /** constraints on the input and output data*/
        std::map<const AbstractData*, std::vector<const Property*> > m_constraints;
    };

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_ABSTRACT_SERVICE_H
/*----------------------------------------------------------------------------*/
