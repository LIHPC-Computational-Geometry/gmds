/*----------------------------------------------------------------------------*/
/** \file    Log.h
 *  \author  F. LEDOUX
 *  \date    09/17/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_LOG_H_
#define GMDS_LOG_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
/*---------------------------------------------------------------------------*/
#include "LogStream.h"
#include "CommonTypes.h"
#include "GMDSUtils_export.h"
/*---------------------------------------------------------------------------*/
namespace gmds {
    /*-----------------------------------------------------------------------*/
    /** \class Log
     *
     *  \brief Class defining a log manager
     */
    /*-----------------------------------------------------------------------*/
    class GMDSUtils_API Log {
        
    public:
        
        static Log& mng();
        static LogLevel& reportingLevel();
        
        virtual ~Log();
        void addStream(LogStream& AStream);
        
        void clear(){out_.clear();}
        /*------------------------------------------------------------------*/
#ifdef DEBUG_GMDS
        /*------------------------------------------------------------------*/
        template<class T> friend inline const Log&
        operator<<(const Log& o,const T& AT)
        {
            for(unsigned int i=0; i<o.out_.size();i++)
            {
                if(o.out_[i].level()<=Log::reportingLevel())
                    o.out_[i] << AT;
            }
            return o;
        }
        
        void flush(){
            for(unsigned int i=0; i<out_.size();i++)
            {
                out_[i].flush();
            }
        }
        /*------------------------------------------------------------------*/
#else
        /*------------------------------------------------------------------*/
        template<class T> friend inline const Log&
        operator<<(const Log& o,const T& AT)
        {
            for(unsigned int i=0; i<o.out_.size();i++)
            {
                if(o.out_[i].level()<=Log::reportingLevel())
                    o.out_[i] << AT;
            }
            return o;

//            //do nothing
//            return o;
        }
        
        void flush(){
            for(unsigned int i=0; i<out_.size();i++)
            {
                out_[i].flush();
            }
        }
        /*------------------------------------------------------------------*/
#endif //DEBUG_GMDS
        /*------------------------------------------------------------------*/
        
    protected:
        
        Log();
        Log (const Log& ALog);
        Log& operator = (const Log&);
        
    private:
        
        static Log* mng_;
        std::vector<LogStream> out_;
        static LogLevel reporting_level_;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif  /* GMDS_LOG_H_ */
/*----------------------------------------------------------------------------*/
