
/*----------------------------------------------------------------------------*/
/** \file    Marks32.h
 *  \author  F. LEDOUX
 *  \date    02/03/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MARKS_32_H_
#define GMDS_MARKS_32_H_
/*----------------------------------------------------------------------------*/
#include "CommonTypes.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class Marks32
 *  \brief Defines boolean marks analoguous to a std::bitstd::set<32> bu using only
 *  	   32 bits whatever the target architecture is.
 */
/*----------------------------------------------------------------------------*/
class Marks32 {
public:

	static const int32_t  MASK_TRUE  = 0xFFFFFFFF;
	static const int32_t  MASK_FALSE = 0x00000000;

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     */
    Marks32(int32_t AMask = MASK_FALSE):mask_(AMask){}

    /*------------------------------------------------------------------------*/
    /** \brief  Copy Constructor.
     *
     *  \param AMarks another 32bits marks.
     */
    Marks32(const Marks32& AMarks):mask_(AMarks.mask_){}

    /*------------------------------------------------------------------------*/
    /** \brief  Operator ==
     *
     *  \param AMarks another boolean marks
     *  \return true if this and AMarks are equal.
     */
    bool operator==(Marks32 AMarks) const { return mask_==AMarks.mask_;}

    /*------------------------------------------------------------------------*/
    /** \brief  Reset all the bits to false
     */
    void reset(){mask_ = MASK_FALSE;}

    /*------------------------------------------------------------------------*/
    /** \brief  Operator !=
     *
     *  \param AMarks another boolean marks
     *  \return true if this and AMarks are different.
     */
    bool operator!=(Marks32 AMarks) const { return mask_!=AMarks.mask_;}

    /*------------------------------------------------------------------------*/
    /** \brief  Operator[]
     *
     *  \param AIndex the index of the bit we want the value
     *  \return the value of AIndex bit.
     */
    bool operator[](const TInt& AIndex) const {return (mask_|(1<<AIndex))==mask_;}

    /*------------------------------------------------------------------------*/
    /** \brief  Sets APos bit to true.
     *
     *  \param APos position of the bit we want to set to true
     */
    void set(int APos){mask_=mask_|(1<<APos);}

    /*------------------------------------------------------------------------*/
    /** \brief  Sets APos bit to false.
     *
     *  \param APos position of the bit we want to set to false
     */
    void unset(int APos){mask_=mask_&(MASK_TRUE^(1<<APos));}

    /*------------------------------------------------------------------------*/
    /** \brief  Sets APos bit to AValue.
     *
     *  \param APos position of the bit we want to set to AValue.
     *  \param AValue the new value of the bit
     */
    void set(int APos,bool AValue) {(AValue)?set(APos):unset(APos);}

    /*------------------------------------------------------------------------*/
    /** \brief  Gets the value of the  APos bit.
     *
     *  \return a boolean indicating the value of APos bit
     */
    bool value(int APos) const {return (mask_|(1<<APos))==mask_;}

    /*------------------------------------------------------------------------*/
    /** \brief  Performs a xor(^) with AMarks
     *
     *  \param AMarks another boolean marks
     *  \return the result of the xor between this and AMarks.
     */
    Marks32 operator^(Marks32 AMarks) const
    { return Marks32(mask_^AMarks.mask_);}

    /*------------------------------------------------------------------------*/
    /** \brief  Inverts the value of the  APos bit.
     *
     *  \param APos the position of the bit we want to invert the value.
     */
    void flip(int pos){mask_=mask_^(1<<pos);}

    /*------------------------------------------------------------------------*/
    /** \brief  Display the all the bit values
     *
     */
    void display() const {
      for(int i=31;i>=0;i--){
    	  (value(i))?std::cout<<"1 ":std::cout<<"0 ";
      }
	  std::cout<<std::endl;
	}

protected:
	/** the mask that store the bits*/
	int32_t mask_;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MARKS_32_H_ */
/*----------------------------------------------------------------------------*/
