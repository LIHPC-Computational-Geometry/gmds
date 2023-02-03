
/*----------------------------------------------------------------------------*/
/** \file    LogStream.h
 *  \author  F. LEDOUX
 *  \date    09/17/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_LOG_STREAM_H_
#define GMDS_LOG_STREAM_H_
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
namespace gmds {
/*---------------------------------------------------------------------------*/
typedef enum {
    LOG_ERROR,  // 0
    LOG_WARNING,// 1
    LOG_INFO,   // 2
    LOG_DEBUG,  // 3
    LOG_DEBUG1, // 4
    LOG_DEBUG2  // 5
} LogLevel;
/*---------------------------------------------------------------------------*/
    class Log;
    /*---------------------------------------------------------------------------*/
    /** \class LogStream
     *
     *  \brief Class defining an output to be logged
     */
    /*----------------------------------------------------------------------------*/
    class LogStream{
        
    public:
        
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor.
         *
         *  \param AStream the stream the LogStream is connected to
         *  \param ALevel  the level of this output
         */
        LogStream(std::basic_ostream<char> *AStream, const LogLevel& ALevel =LOG_INFO)
        :m_stream(AStream),m_level(ALevel){}
        
        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor.
         *
         *  \param ALog the LogStream to copy
         */
        LogStream(const LogStream& ALog) {
            m_stream = ALog.m_stream;
            m_level=ALog.m_level;
        }
        
        
        ~LogStream(){}
        
        /*------------------------------------------------------------------------*/
        /** \brief Overloading of =
         *
         *  \param ALog the LogStream to copy
         */
        const LogStream& operator=(const LogStream& ALog){
            m_stream = ALog.m_stream;
            m_level  = ALog.m_level;
            return *this;
        }
        
        /*------------------------------------------------------------------------*/
        /** \brief give access to the stream
         *
         *  \return the stream
         */
        inline std::ostream& stream() const {return *m_stream;}
        
        /*------------------------------------------------------------------------*/
        /** \brief give consultation access to the level
         *
         *  \return the stream
         */
        inline LogLevel   level() const {return m_level;}
        /*------------------------------------------------------------------------*/
        /** \brief give modification access to the level
         *
         *  \return the stream
         */
        inline LogLevel& level() {return m_level;}
        
        /*----------------------------------------------------------------------------*/
        template<class T> friend inline const LogStream&
        operator<<(const LogStream& ALO,  const T& AT)
        { ALO.stream() << AT; return ALO;}
        
        void flush(){m_stream->flush();}
    protected:
        
        std::ostream* m_stream;
        LogLevel      m_level;
        
    };
    
    /*---------------------------------------------------------------------------*/
    /** \class FileLogStream
     *
     *  \brief Class defining an file stream to be logged
     */
    /*----------------------------------------------------------------------------*/
    class FileLogStream: public LogStream{
        
    public:
        
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor.
         *
         *  \param AFileName  the name of the file
         *  \param ALevel  the level of this output
         */
        FileLogStream(const std::string& AFileName,
                      const LogLevel& ALevel =LOG_INFO)
        :LogStream(NULL,ALevel) {
            m_file_name = AFileName;
            m_file_stream = new std::ofstream(m_file_name.c_str(),
                                              std::ios::out | std::ios::app);
            m_stream = m_file_stream;
        }
        
        
        
        ~FileLogStream() {
            m_file_stream->close();
            if(m_file_stream!=NULL)
                delete  m_file_stream;
        }
        
    private:
        std::ofstream* m_file_stream;
        std::string m_file_name;
        
        
    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif  /* GMDS_LOG_STREAM_H_ */
/*----------------------------------------------------------------------------*/
