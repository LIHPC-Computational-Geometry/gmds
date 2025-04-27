/*---------------------------------------------------------------------------*/
#ifndef MATCHING_IACTION_H
#define MATCHING_IACTION_H
/*---------------------------------------------------------------------------*/
#include <memory>
/*---------------------------------------------------------------------------*/
class IState;
/*---------------------------------------------------------------------------*/
struct IAction {

    /**@brief Computes the state reach from @p AState by applying the current
     * action
     * @param[in] AState the state we start from
     * @return the state that is built from this @p AState when applying the
     * action
     */
    virtual std::shared_ptr<IState> apply_on(std::shared_ptr<IState> AState) const =0;

    virtual ~IAction() = default;
    virtual bool operator==(const IAction& other) const = 0;
	 virtual std::string get_description() const = 0;
};
/*---------------------------------------------------------------------------*/
#endif //MATCHING_IACTION_H
/*---------------------------------------------------------------------------*/
