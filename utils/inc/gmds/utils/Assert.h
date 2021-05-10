/*----------------------------------------------------------------------------*/
/** \file    Assert.h
 *  \author  F. LEDOUX
 *  \date    04/11/2016
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ASSERT_H_
#define GMDS_ASSERT_H_
/*----------------------------------------------------------------------------*/
#include "CommonFlags.h"
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    /**
     * \enum  Assert termination mode
     * \brief Defines how assertions must be terminated, either with launching
     *        an exception or abording the program.
     */
    enum AssertionMode {
        ASSERT_MODE_THROW,
        ASSERT_MODE_ABORT
    };
    /*------------------------------------------------------------------------*/
    /** \brief Sets the assertion mode to \p AMode.
     *
     *  \param[in] AMode the assertion mode
     */
    void EXPORT_GMDS setAssertMode(AssertionMode AMode);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the current assertion mode
     */
    AssertionMode EXPORT_GMDS getAssertMode();

    /*------------------------------------------------------------------------*/
    /**
     * \brief Prints an assertion failure when \p ACondition is evaluated to 
     *        false. The termination is done according to the current assertion
     *        mode.
     *
     * \param[in] ACondition String representation of the condition
     * \param[in] AFileName  Name of the file where the assertion is raised
     * \param[in] ALineNum   Number Line where the assertion is raised
     */
    void EXPORT_GMDS assertFailed(const std::string& ACondition,
                                  const std::string& AFileName,
                                  int ALineNum);
    
    /*------------------------------------------------------------------------*/
    /** \brief Prints a range assertion failure when \p AVal is out of the range
     *         [\p AMin, \p AMax\]. The termination is done according to the
     *         current assertion mode.
     *
     * \param[in] AVal      the value to be checked
     * \param[in] AMin      Minimum value
     * \param[in] AMax      Maximum value
     * \param[in] AFileName Name of the file where the assertion is raised
     * \param[in] ALineNum  Number Line where the assertion is raised     
     */
    void EXPORT_GMDS assertRangeFailed(const double& AVal,
                                       const double& AMin,
                                       const double& AMax,
                                       const std::string& AFileName,
                                       int ALineNum);
}
/*----------------------------------------------------------------------------*/
// Two levels of assert,
// use gmds_XXX_assert()        for assertions called in release and debug mode
// use gmds_XXX__debug_assert() for assertions called only in debug mode
/*----------------------------------------------------------------------------*/
/** \brief Check that condition \p AX is verified (so evaluated to true)
 * 
 * \param[in] AX the boolean expression of the condition
 *
 * \see gmds::assert_failed()
 */
#define GMDS_ASSERT(AX) {                                \
        if(!(AX)) {                                      \
            gmds::assertFailed(#AX, __FILE__, __LINE__);\
        }                                                \
}
/*----------------------------------------------------------------------------*/
/** \brief Check that that \p AV is in [\p AMin, \p AMax]
 *
 * \param[in] AX   the value to be checked
 * \param[in] AMin minimum authorized value
 * \param[in] AMax maximum authorized value
 *
 * \see gmds::assert_range_failed()
 */
#define GMDS_RANGE_ASSERT(AX, AMin, AMax) {        \
        if(((AX) < (AMin)) || ((AX) > (AMax))) {   \
            gmds::assertRangeFailed(AX, AMin, AMax \
                __FILE__, __LINE__                 \
            );                                     \
        }                                          \
}
/*----------------------------------------------------------------------------*/
#endif  /* GMDS_ASSERT_H_ */
/*----------------------------------------------------------------------------*/


