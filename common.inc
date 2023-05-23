// common.inc

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Common macros
*   To be included into .cpp-files
*
*/


#ifndef COMMON_INC
#define COMMON_INC


#define FOR_START(type,i,start,end)  for (type i = (start); i < (end); i++)  
#define FOR(type,i,end)  FOR_START (type, i, 0, (end))
// Fast FOR: precomputed end
#define FFOR_START(type,i,start,end)  for (type i = (start), i ## _end_ = (end); i < i ## _end_; i++)  
#define FFOR(type,i,end)  FFOR_START (type, i, 0, (end))
//
#define FOR_REV_END(type,i,end,start)  for (type i = (start); i-- > (end);)  
#define FOR_REV(type,i,start)  FOR_REV_END(type, i, 0, (start))
// FOR(type,i,a) and FOR_REV(type,i,a) iterate over the same i
// FOR_START(type,i,a,b) and FOR_REV_END(type,i,a,b) iterate over the same i

#define ITER(ContainerType,iter,container) \
  for (ContainerType::iterator iter = (container). begin (); iter != (container). end (); iter++)
#define CONST_ITER(ContainerType,iter,container) \
  for (ContainerType::const_iterator iter = (container). begin (); iter != (container). end (); iter++)
#define ITER_REV(ContainerType,iter,container) \
  for (ContainerType::reverse_iterator iter = (container). rbegin (); iter != (container). rend (); iter++)
#define CONST_ITER_REV(ContainerType,iter,container) \
  for (ContainerType::const_reverse_iterator iter = (container). rbegin (); iter != (container). rend (); iter++)


#define  Case     break; case
#define  Default  break; default


// Exceptions

#include <string>
#include <exception>


#ifdef _MSC_VER
  #define FUNC  std::string (__FUNCSIG__) + ":\n"
#else
  #define FUNC  std::string (__PRETTY_FUNCTION__) + ":\n"
#endif


namespace Common_sp
{
  [[noreturn]] void throwf (const std::string &s);
}

#define ERROR_MSG(msg)                     \
  { if (! std::uncaught_exceptions ())     \
      Common_sp::throwf (std::string ("\"" __FILE__ "\", line ") + to_string (__LINE__) + ": " + (msg));  \
    else {}                                \
    exit (1);                              \
  }
#define ERROR            ERROR_MSG (FUNC "ERROR")
#define NOT_IMPLEMENTED  ERROR_MSG (FUNC "NOT IMPLEMENTED")
#define NEVER_CALL       ERROR_MSG (FUNC "NEVER CALL")


#define QC_ASSERT(cond)  if (! (cond)) ERROR_MSG (#cond) else {}   


// Logic errors
#ifdef NDEBUG                                       
  #define ASSERT(cond)
#else
  #define ASSERT(cond)  QC_ASSERT (cond)
#endif


#ifdef NDEBUG                                       
  #define EXEC_ASSERT(cond)  cond
#else
  #define EXEC_ASSERT(cond)  { const bool c_ = (cond); if (! c_) ERROR_MSG (#cond); }
#endif


#define IMPLY(a,b)     { if (a) { ASSERT (b) }}
#define QC_IMPLY(a,b)  { if (a) { QC_ASSERT (b) }}

#define ASSERT_EQ(x,y,delta)     ASSERT    (std::fabs((x) - (y)) <= (delta))
#define QC_ASSERT_EQ(x,y,delta)  QC_ASSERT (std::fabs((x) - (y)) <= (delta))



#define MODULE_INIT   static bool run = false; \
                      if (run)                 \
                        return true;           \
                      run = true;


#define LESS_PART(x,y,part)  { if ((x).part < (y).part) return true; \
                               if ((y).part < (x).part) return false; }
#define LESS_COMP(comp)  { const ebool c = comp; if (c == etrue) return true; if (c == efalse) return false; }



#define PRINT(x)  { Offset::newLn (cout); cout << #x << " = " << (x); }


#define XSTR(s) STR(s)
#define STR(s)  #s



#endif
