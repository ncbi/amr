// common.hpp

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
*   Common utilities
*
*/


#ifndef COMMON_HPP_64052  // random number
#define COMMON_HPP_64052


#ifdef _MSC_VER
  #pragma warning (disable : 4061)  // enumerator ... in switch of enum ... is not explicitly handled by a case label
  #pragma warning (disable : 4290)  // C++ exception specification ignored except to indicate a function is not __declspec(nothrow)
  #pragma warning (disable : 4365)  // conversion from 'type_1' to 'type_2', signed/unsigned mismatch (bool -> size_t)
  #pragma warning (disable : 4503)  // decorated name length exceeded, name was truncated
  #pragma warning (disable : 4514)  // '...': unreferenced inline function has been removed
  #pragma warning (disable : 4521)  // multiple copy constructors specified
  #pragma warning (disable : 4522)  // multiple assignment operators specified
  #pragma warning (disable : 4592)  // symbol will be dynamically initialized (implementation limitation)
  #pragma warning (disable : 4625)  // copy constructor was implicitly defined as deleted
  #pragma warning (disable : 4626)  // assignment operator was implicitly defined as deleted
  #pragma warning (disable : 4710)  // function not inlined
  #pragma warning (disable : 4800)  // 'const char *': forcing value to bool 'true' or 'false' (performance warning)
  #pragma warning (disable : 4820)  // '...' bytes padding added after data member '...'
  #pragma warning (disable : 5026)  // move constructor was implicitly defined as deleted
  #pragma warning (disable : 5027)  // move assignment operator was implicitly defined as deleted

  #pragma warning (disable : 4005)  // macro redefinition
  #define _HAS_ITERATOR_DEBUGGING 0
  #pragma warning (default : 4005)  
#endif

#include <cassert>
#include <cstring>
#include <cmath>
#include <string>
#include <stdexcept>
#include <limits>
#include <array>
#include <list>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <filesystem>

#include <thread>
#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable:4265)
#endif
#include <mutex>
#ifdef _MSC_VER
	#pragma warning(pop)
#endif

using namespace std;



namespace Common_sp
{


bool initCommon ();
  // Module initialization
  // Invoked automaticallly


// Numeric types

typedef  unsigned char  uchar; 
typedef  unsigned int   uint; 
typedef  unsigned long  ulong; 



// Numeric constants

constexpr size_t no_index = numeric_limits<size_t>::max ();
static_assert ((size_t) 0 - 1 == no_index);

constexpr double NaN = numeric_limits<double>::quiet_NaN ();  
  


// Global variables

extern vector<string> programArgs;

extern string programName;

string getCommandLine ();

extern ostream* logPtr;

extern bool qc_on;

extern ulong seed_global;
  // >= 1



// Errors

constexpr const char* error_caption ("*** ERROR ***");

[[noreturn]] void errorExit (const char* msg,
                             bool segmFault = false);
  // Input: programArgs, programName, logptr
	// Update: *logPtr
	// Invokes: if segmFault then abort() else exit(1)

void errorExitStr (const string &msg);
  // For debugger: should not be inline
  // Invokes: beep()

[[noreturn]] void throwf (const string &s); 
  // For debugger: should not be inline
  // Invokes: throw logic_error

inline void reportException (const string &msg,
                             const exception &e,
                             bool force)
  { const string s (msg + ": " + e. what ());
    if (force)
      cerr << s << endl;
    else
      throw runtime_error (s);
  }    	    	


void sleepNano (long nanoSec);



void beep ();
  // Requires: !isRedirected()
  //           SHLVL = 1 ??
    


template <typename T> 
  inline bool isNan (T x)
    { return x != x; }


// Comparison templates

template <typename T> 
  inline void swapGreater (T &a, T &b)
    { if (a > b)
        swap (a, b);
    }

template <typename T> 
  inline bool maximize (T &a, T b)
    { if (a < b) { a = b; return true; } return false; }

template <typename T> 
  inline bool minimize (T &a, T b)
    { if (a > b) { a = b; return true; } return false; }
    	
template <typename T>
  inline T difference (T a, T b)
    { if (a > b) return a - b; return b - a; }

template <typename T /*:number*/> 
  inline bool between (T x, T low, T high)
    { return x >= low && x < high; }

template <typename T/*:number*/> 
  inline bool betweenEqual (T x, T low, T high)
    { return x >= low && x <= high; }

template <typename T>
  inline bool lessPtr (const T* x,
                       const T* y)
    { return *x < *y; }



typedef int (*CompareInt) (const void*, const void*);



// Pointer templates

template <typename T>
  inline T* var_cast (const T* t)
    { return const_cast <T*> (t); }
template <typename T>
  inline T& var_cast (const T &t)
    { return const_cast <T&> (t); }

template <typename T, typename S>
  inline T const_static_cast (const S* s)
    { return static_cast <T> (const_cast <S*> (s)); }
template <typename T, typename S>
  inline T const_static_cast (const S &s)
    { return static_cast <T> (const_cast <S&> (s)); }

template <typename T>
  inline T* checkPtr (T* t)
    { if (t)
        return t;
      throwf ("Dereferencing nullptr");
      return nullptr;
    }

  

// Container templates

template <typename T>
  inline void sort (T &t)
    { std::sort (t. begin (), t. end ()); }

template <typename T, typename StrictlyLess>
  inline void sort (T &t,
                    const StrictlyLess &strictlyLess)
    { std::sort (t. begin (), t. end (), strictlyLess); }

template <typename T, typename U>
  bool intersects (const T &t,
                   const U &u)
    { if (u. empty ())
        return false;
      for (const auto& x : t)
        if (u. find (static_cast <typename U::value_type> (x)) != u. end ())
          return true;
      return false;
    }

template <typename T>
  unordered_set<T> setIntersect (const unordered_set<T> &s1,
                                 const unordered_set<T> &s2)
    { const unordered_set<T>* s1_ = & s1;
      const unordered_set<T>* s2_ = & s2;
      if (s1. size () > s2. size ())
        swap (s1_, s2_);
      unordered_set<T> res;  res. rehash (s1_->size ());
      for (const T& t : *s1_)
        if (contains (*s2_, t))
          res. insert (t);
      return res;
    }
    
template <typename T, typename U>
  bool containsSubset (const T &t,
                       const U &u)
    { for (const auto& x : u)
        if (t. find (static_cast <typename U::value_type> (x)) == t. end ())
          return false;
      return true;
    }

template <typename T>
  unordered_set<T> getSetMinus (const unordered_set<T> &s1,
                                const unordered_set<T> &s2)
    { unordered_set<T> res;  res. rehash (s1. size ());
      for (const T& t : s1)
        if (! contains (s2, t))
          res. insert (t);
      return res;
    }

template <typename T, typename U>
  void setMinus (T &t,
                 const U &u)
    { for (const auto& x : u)
        t. erase (x);
    }

template <typename Key, typename Value, typename KeyParent>
  inline bool contains (const map <Key, Value> &m,
                        const KeyParent& key)
    { return m. find (key) != m. end (); }

template <typename Key, typename Value, typename KeyParent>
  inline bool contains (const unordered_map <Key, Value> &m,
                        const KeyParent& key)
    { return m. find (key) != m. end (); }

template <typename Key, typename KeyParent>
  inline bool contains (const unordered_set <Key> &m,
                        const KeyParent& key)
    { return m. find (key) != m. end (); }

template <typename T, size_t N>
  size_t indexOf (const array<T,N> &arr, const T item)
    { for (size_t i = 0; i < N; i++)
        if (arr [i] == item)
          return i;
      return no_index;
    }

template <typename T, size_t N>
  inline bool contains (const array<T,N> &arr, const T item)
    { return indexOf (arr, item) != no_index; }

template <typename Key, typename Value, typename KeyParent>
  bool find (const map <Key, Value> &m,
             const KeyParent& key,
             Value &value)
    // Return: success
    // Output: value, if Return
    { const auto& it = m. find (key);
    	if (it == m. end ())
    		return false;
    	value = it->second; 
    	return true;
    }

template <typename Key, typename Value, typename KeyParent>
  bool find (const unordered_map <Key, Value> &m,
             const KeyParent& key,
             Value &value)
    // Return: success
    // Output: value, if Return
    { const auto& it = m. find (key);
    	if (it == m. end ())
    		return false;
    	value = it->second; 
    	return true;
    }

template <typename Key, typename Value, typename Hash, typename KeyParent>
  bool find (const unordered_map <Key, Value, Hash> &m,
             const KeyParent& key,
             Value &value)
    // Return: success
    // Output: value, if Return
    { const auto& it = m. find (key);
    	if (it == m. end ())
    		return false;
    	value = it->second; 
    	return true;
    }

template <typename Key, typename Value, typename KeyParent>
  const Value* findPtr (const map <Key, Value> &m,
                        const KeyParent& key)
    { const auto& it = m. find (key);
    	if (it == m. end ())
    		return nullptr;
    	return & it->second; 
    }

template <typename Key, typename Value, typename KeyParent>
  const Value* findPtr (const unordered_map <Key, Value> &m,
                        const KeyParent& key)
    { const auto& it = m. find (key);
    	if (it == m. end ())
    		return nullptr;
    	return & it->second; 
    }

template <typename Key, typename Value, typename KeyParent>
  const Value* findPtr (const map <Key, const Value* /*!nullptr*/> &m,
                        const KeyParent& key)
    { const Value* value;
      if (find (m, key, value))
        return value;
      return nullptr;
    }

template <typename Key, typename Value, typename KeyParent>
  const Value* findPtr (const unordered_map <Key, const Value* /*!nullptr*/> &m,
                        const KeyParent& key)
    { const Value* value;
      if (find (m, key, value))
        return value;
      return nullptr;
    }

template <typename Key, typename Value>
  const Value& findMake (map <Key, const Value* /*!nullptr*/> &m,
                         const Key& key)
    { if (! contains (m, key))
        m [key] = new Value ();
      return * m [key];
    }

template <typename Key, typename Value>
  const Value& findMake (unordered_map <Key, const Value* /*!nullptr*/> &m,
                         const Key& key)
    { if (! contains (m, key))
        m [key] = new Value ();
      return * m [key];
    }

template <typename T, typename UnaryPredicate>
  inline size_t count_if (T &t, UnaryPredicate pred)
    { return (size_t) std::count_if (t. begin (), t. end (), pred); }

template <typename To, typename From>
  inline void insertAll (To &to,
                         const From &from)
    { 
    #pragma GCC diagnostic ignored "-Waggressive-loop-optimizations"
      to. insert (to. begin (), from. begin (), from. end ()); 
    #pragma GCC diagnostic warning "-Waggressive-loop-optimizations"
    }

template <typename To, typename From>
  inline void insertIter (To &to,
                          From &&from)
    { for (auto&& x : from)
        to. insert (std::move (x));
    }

template <typename To, typename From>
  inline void insertIter (To &to,
                          const From &from)
    { for (const auto& x : from)
        to. insert (x);
    }

template <typename From, typename What>
  inline void eraseIter (From &from,
                         const What &what)
    { for (const auto& x : what)
        if (contains (from, x))
          from. erase (x);
    }

template <typename T /*container*/>
  void save (ostream &os,
             const T &container,
             char sep)
    { bool first = true;
      for (const auto& t : container)
      { if (! first)
          os << sep;
        os << t;
        first = false;
      }
    }



// KeyValue

typedef  map<string/*key*/,string/*value*/>  KeyValue;

inline string find (const KeyValue& kv,
                    const string& key,
                    bool force)
  { const auto& it = kv. find (key);
    if (it == kv. end ())
    { if (force)
        return key;
      throw runtime_error ("Key \"" + key + "\" is not found");
    }
    return it->second;
  }



// hash

extern hash<string> str_hash;

extern hash<size_t> size_hash;

// PAR
constexpr size_t small_hash_class_max = 1000;  
constexpr size_t large_hash_class_max = small_hash_class_max * 100;  

inline size_t str2hash_class (const string &s,
                              bool large_hash)
  { return str_hash (s) % (large_hash 
                             ? large_hash_class_max 
                             : small_hash_class_max
                          ); 
  }
 


// bool

inline void toggle (bool &b)
  { b = ! b; }

inline int getSign (bool b)
  { return b ? 1 : -1; }

inline bool boolPow (bool x, bool power)
  { return power ? x : ! x; }
  	
inline string yesNo (bool x)
  { return x ? "Y" : "N"; }



// ebool: extended bool

enum ebool {efalse = false, 
            etrue = true, 
            enull = true + 1};

inline ebool toEbool (bool b)
  { return b ? etrue : efalse; }


inline bool operator<= (ebool a, ebool b)
  { static constexpr char rank [3/*ebool*/] = {0, 2, 1};
  	return rank [a] <= rank [b];
  }

inline ebool negate (ebool b)
  { switch (b)
    { case efalse: return etrue;
      case etrue:  return efalse;
      default:     return b;
    }
  }

inline void toggle (ebool &b)
  { b = negate (b); }
  
inline string ebool2txt (ebool choice,
                         const string &yes,
                         const string &no,
                         const string &ambig)
  { switch (choice)
  	{ case etrue:  return yes;
  		case efalse: return no;
  		default:     return ambig;
    }
  }



// char

inline bool isChar (long long n)
  { return between<long long> (n, -128, 128); }

inline bool isAlpha (char c)
  { return strchr ("abcdefghijklmnopqrstuvwxyz", tolower (c)); }
  // isalpha() is locale-specific

inline bool isDigit (char c)
  { return strchr ("0123456789", c); }
  // isdigit() is locale-specific
  
inline bool isLetter (char c)
  { return isAlpha (c) || isDigit (c) || c == '_'; }

inline char toUpper (char c)
  { return (char) toupper (c); }

inline char toLower (char c)
  { return (char) tolower (c); }

inline bool isUpper (char c)
  { return toUpper (c) == c; }

inline bool isLower (char c)
  { return toLower (c) == c; }

inline bool printable (char c)
  { return between (c, ' ', (char) 127); }
  
inline bool nonPrintable (long c)
  { return between<long> (c, '\x001', ' '); }
  
string nonPrintable2str (char c);

inline bool isDelimiter (char c)
  { return    c != ' '
           && printable (c)
           && ! isLetter (c);
  }
  
inline bool isSpace (char c)
  { return c > '\0' && c <= ' ' && isspace (c); }

string to_url (const string &s);


// char*

inline const char* nvl (const char* s,
                        const char* nullS = "-")
  { return s ? s : nullS; }
  	
  	  	

// string

extern const string noString;

#ifndef _MSC_VER
  inline string getEnv (const string& name)
    { if (const char* p = getenv (name. c_str ()))
        return p;
      return noString;
    }
#endif
  
inline string ifS (bool cond,
                   const string &s)
  { return cond ? s : noString; }

inline string nvl (const string& s,
                   const string& nullS = "-")
  { return s. empty () ? nullS : s; }
  	
inline string appendS (const string &s,
                       const string &suffix)
  { return s. empty () ? noString : (s + suffix); }
  	
inline string prependS (const string &s,
                        const string &prefix)
  { return s. empty () ? noString : (prefix + s); }
  	
inline void add (string &to,
                 const string &delimiter,
		             const string &what)
  { if (! to. empty ())
  	  to += delimiter;
  	to += what;
  }

inline bool isQuoted (const string &s,
                      char quote = '\"')
  { return ! s. empty () && s [0] == quote && s [s. size () - 1] == quote; }

string strQuote (const string &s,
                 char quote = '\"');

inline string unQuote (const string &s)
  { return s. substr (1, s. size () - 2); }

bool strBlank (const string &s);

bool getScientific (string numberS,
                    bool &hasPoint,
                    streamsize &decimals);

template <typename T>
  string toString (const T &t)
    { ostringstream oss;
      oss << t;
      return oss. str ();
    }

template <typename T>
  T str2 (const string &s)
    { static_assert (numeric_limits<T>::max() > 256, "str2 does not work on chars");
    	T i;
    	istringstream iss (s);
      iss >> i;
      if (   Common_sp::strBlank (s)
          || ! iss. eof ()
          || iss. fail ()
         )
        throw runtime_error ("Cannot convert " + strQuote (s) + " to number");
      return i;
    }

template <typename T>
  bool str2 (const string &s,
             T &t)
    { try { t = str2<T> (s); return true; } 
        catch (...) { return false; } 
    }

void commaize (string &s);
  // ' ' --> ','

inline bool isLeftBlank (const string &s,
                         size_t spaces)
  { size_t i = 0;
    while (i < s. size () && i < spaces && s [i] == ' ')
      i++;
    return i == spaces;
  }

string pad (const string &s,
            size_t size,
            ebool right,
            char space = ' ');
  // right: enull - center
  
template <typename T>
  string padNatural (T natural,
                     size_t size)
    { return pad (to_string (natural), size, etrue, '0'); }

bool goodName (const string &name);

bool isIdentifier (const string& name,
                   bool dashInName);
  // Return: true. !empty(), !dashInName => (c \in mame => isLetter(c))

bool isNatural (const string& name,
                bool leadingZeroAllowed);

void strUpper (string &s);

void strLower (string &s);

bool isUpper (const string &s);

bool isLower (const string &s);

inline string strUpper1 (const string &s)
  { if (s. empty ())
      return s;
    return toUpper (s [0]) + s. substr (1);
  }

inline bool charInSet (char c,
		                   const string &charSet)
  { return charSet. find (c) != string::npos; }

size_t stringNotInSet (const string &s,
                   	   const string &charSet,
                   	   ebool uppercase);
  // Return: index in s which is not in charSet; may be no_index

size_t strCountSet (const string &s,
		                const string &charSet);

void strDeleteSet (string &s,
		               const string &charSet);

inline string getLeft (const string &s,
                       size_t size)
  { return s. substr (0, size); }

inline string getRight (const string &s,
                        size_t size)
  { if (size >= s. size ())
      return s;
    return s. substr (s. size () - size);
  }

inline bool isLeft (const string &s,
                    const string &left)
  { return getLeft (s, left. size ()) /*s. substr (0, left. size ())*/ == left; }

inline bool isRight (const string &s,
                     const string &right)
  { return getRight (s, right. size ()) == right; }

bool trimPrefix (string &s,
                 const string &prefix);
  // Return: success

bool trimSuffix (string &s,
                 const string &suffix);
  // Return: success

void trimSuffixNonAlphaNum (string &s);

bool trimTailAt (string &s,
                 const string &tailStart);
  // Return: trimmed; tailStart.empty() => false
  // Update: s

void trimLeading (string &s);

void trimTrailing (string &s);

inline void trim (string &s)
  { trimTrailing (s);
    trimLeading (s); 
  }

void trimLeading (string &s,
                  char c);
  // Invokes: isSpace()

void trimTrailing (string &s,
                   char c);
  // Invokes: isSpace()

inline void trim (string &s,
                  char c)
  { trimTrailing (s, c);
    trimLeading  (s, c); 
  }

inline bool strNull (const string &s)
  { if (strBlank (s))
      return true;
    string s1 (s);
    trim (s1);
    strUpper (s1);
    return    s1 == "NULL"
           || s1 == "NA"
           || s1 == "N/A";
  }

inline bool contains (const string &hay,
                      const string &needle)
  { return hay. find (needle) != string::npos; }

inline bool contains (const string &hay,
                      char needle)
  { return hay. find (needle) != string::npos; }

size_t containsWord (const string& hay,
                     const string& needle);
  // Return: position of needle in hay; string::npos if needle does not exist
  
struct StringMatch
{
  enum Type {part, word, whole};
};

bool matches (const string &hay,
              const string &needle,
				      StringMatch::Type type);

void replace (string &s,
              char from,
              char to);
  
void replace (string &s,
              const string &fromAnyChars,
              char to);
  
void replaceStr (string &s,
                 const string &from,
                 const string &to);
  // Replaces <from> by <to> in s from left to right
  // If <to> contains <from> then the replacing <to> is skipped else postcondtion is that s does not contain <from>
  // Requires: !from.empty()

string to_c (const string &s);
  // " --> \", etc.

void collapseSpace (string &s);
  // "  " --> " "
  
void visualizeTrailingSpaces (string &s);
  // Invokes: nonPrintable2str(' ')
  
string str2streamWord (const string &s,
                       size_t wordNum);
  // Return: May be empty()
  // Input: wordNum: 0-based  

string str2sql (const string &s);
  // Return: ' s ' (replacing ' by '')

string findSplit (string &s,
                  char c = ' ');
	// Return: prefix of s+c before c
	// Update: s

string rfindSplit (string &s,
                   char c = ' ');
	// Return: suffix of c+s after c
	// Update: s

void reverse (string &s);

size_t strMonth2num (const string& month);
  // Input: month: "Jan", "Feb", ... (3 characters)



// Bit operations

inline bool contains (uint who,
                      uint what)
  { return (who & what) == what; }

inline uchar reverse (uchar b)
  { return uchar (numeric_limits<uchar>::max () - b); }

size_t byte2first (uchar b);
  // Return: First 1-bit, 0-based

size_t utf8_len (char first);
  // Return: 0 - first is ASCII
  //         1 - first is UTF-8 not first byte
  //         else - first is UTF-8 first byte and value is the number of UTF-8 bytes



// Integer operations

inline void advance (size_t &index, 
                     size_t size)
  // Update: index: < size
  { index++; if (index == size) index = 0; }
  
template <typename T/*:integer*/> 
  inline bool even (T x)
    { static_assert (numeric_limits<T>::is_integer);
      return x % 2 == 0; 
    }

inline bool divisible (size_t n,
                       size_t divisor)
  { return ! (n % divisor); }
  
inline uint remainder (int n, uint div)
  { if (const int rem = n % (int) div) 
      return (uint) (rem > 0 ? rem : (rem + (int) div)); 
    return 0; 
  }
  
inline uint gcd (uint a,
                 uint b)
	// Greatest common divisor
	// Euclid's algorithm
	{ if (! b)
			return a;
		if (a < b)
			return gcd (b, a);
		return gcd (b, a % b);
	}

size_t powInt (size_t a,
               size_t b);
  // Return: a^b
  // Time: O(log(b))



// Time

string getNow ();
  // Return: YYYY-MM-DD HH:mm:SS



struct Rand
// Numerical Recipes in C, p. 279
{
	static const long /*actually ulong*/ max_;
private:
	long seed {0};
	  // 0 < seed < max_
public:
	

	explicit Rand (ulong seed_arg = 1)
  	{ setSeed (seed_arg); }
	void setSeed (ulong seed_arg)
    {	seed = (long) (seed_arg % (ulong) max_);
	    qc ();
    }
    // Input: seed_arg > 0

	  
	ulong get (ulong max);
    // Return: in 0 .. max - 1
    // Input: 0 < max <= max_
	ulong get ()
	  { return get ((ulong) max_); }
  double getProb ()
    { return (double) get ((ulong) max_) / ((double) max_ - 1); }
    // Return: [0,1]
private:
	void qc () const;
	void run ();	
};



// Hexadecimal

inline bool isHex (char c)
  { return isDigit (c) || strchr ("ABCDEF", toUpper (c)); }

inline uchar hex2uchar (char c)
  { return uchar (isDigit (c) ? (c - '0') : (toUpper (c) - 'A' + 10)); }

inline string uchar2hex (uchar c)
  { constexpr char hex [16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
    string res ("  ");
    res [0] = hex [c / 16];
    res [1] = hex [c % 16];
    return res;
  }
 
string unpercent (const string &s);
  // '%HH' -> char



struct DisjointCluster
// Cormen, Leiserson, Rivest, Introduction to Algorithms, p. 449
{
protected:
  DisjointCluster* parentDC {nullptr};
    // !nullptr
    // Tree
    // = this <=> root
  size_t rankDC {0};
    // Upper bound on the height of *this
    // (Height = max. # arcs between *this and a leaf)
public:

protected:
	DisjointCluster ()
	  { init (); }
public:
	
	void init ()
    { rankDC = 0;
			parentDC = this;
		}
	void merge (DisjointCluster &other);
  DisjointCluster* getDisjointCluster ();
};



// Simple classes

class Notype {};



struct Nocopy
{
protected:
	Nocopy () = default;
  explicit Nocopy (const Nocopy &) = delete;
  explicit Nocopy (Nocopy &&) = delete;
  Nocopy& operator= (const Nocopy &) = delete;
};



template <typename T>
  struct Singleton : Nocopy
  {
  private:
  	static bool beingRun;
  protected:
  	Singleton ()
  	  { if (beingRun)
  	  	  throwf ("Singleton");
  	  	beingRun = true;
  	  }
   ~Singleton () 
      { beingRun = false; }
  };
template <typename T> bool Singleton<T>::beingRun = false;



template <typename T>
  class Keep : Nocopy
  {
    T* ptr;
    const T t;
    
  public:
    explicit Keep (T &t_arg)
      : ptr (& t_arg)
      , t (t_arg)
      {}
   ~Keep () 
      { if (ptr)
      	  *ptr = t; 
      }
      
    const T& get () const
      { return t; }
    void reject ()
      { ptr = nullptr; }
  };



template <typename T>
  struct Pair : pair <T, T>
  {
  private:
  	typedef  pair <T, T>  P;
  public:
    
    Pair (const T &a,
          const T &b)
      : P (a, b)
      {}
    Pair (T &&a,
          T &&b)
      : P (std::move (a), std::move (b))
      {}
    Pair () = default;
    Pair (Pair<T> &&other)
      : P (std::move (other. first), std::move (other. second))
      {}
    Pair (const Pair<T> &other) = default;
    Pair& operator= (const Pair<T> &other) = default;

      
    bool same () const
      { return P::first == P::second; }
    bool has (const T &t) const
      { return    P::first  == t 
               || P::second == t;
      }
    T& findOther (const T &t) const
      { if (P::first == t)
          return P::second;
        return P::first;
      }
    void swap ()
      { std::swap (P::first, P::second); }
    bool empty () const
      { return    P::first  == T ()
               && P::second == T ();
      }
  };



template <typename T, typename U>
  inline ostream& operator<< (ostream &os,
                              const pair<T,U> &p) 
    { os << p. first << '\t' << p. second;
      return os;
    }



template <typename T>
  struct List : list<T>
  {
  private:
  	typedef  list<T>  P;
  public:


  	List () = default;
    template <typename U/*:<T>*/>
      List (const list<U> &other)
        { *this << other; }
    template <typename U/*:<T>*/>
      explicit List (const vector<U> &other)
        { *this << other; }
    bool empty () const
      {
      #pragma GCC diagnostic ignored "-Wnull-dereference"
        return P::empty ();
      #pragma GCC diagnostic warning "-Wnull-dereference"
      }
    size_t size () const
      {
      #pragma GCC diagnostic ignored "-Wnull-dereference"
        return P::size ();
      #pragma GCC diagnostic warning "-Wnull-dereference"
      }

  	  
  	T at (size_t index) const
  	  { size_t i = 0;
  	  	for (const T& t : *this)
  	  	  if (i == index)
  	  	  	return t;
  	  	  else
  	  	  	i++;
  	  	throw runtime_error ("List index is out of range");
  	  	return T ();  // dummy
  	  }
  	size_t find (const T &t) const
  	  { size_t i = 0;
  	  	for (const T& item : *this)
  	  	  if (item == t)
  	  	  	return i;
  	  	  else
  	  	  	i++;
  	  	return no_index;
  	  }
  	size_t find (const T* t) const
  	  { size_t i = 0;
  	  	for (const T& item : *this)
  	  	  if (& item == t)
  	  	  	return i;
  	  	  else
  	  	  	i++;
  	  	return no_index;
  	  }
    bool isPrefix (const List<T> &prefix) const
      { typename List<T>::const_iterator wholeIt  =      P::begin ();
      	typename List<T>::const_iterator prefixIt = prefix. begin ();
      	for (;;)
      	{ if (prefixIt == prefix. end ())
  	    		return true;
  	    	if (wholeIt == P::end ())
  	    		return false;
  	      if (*prefixIt != *wholeIt)
  	      	return false;
  	      wholeIt ++;
  	      prefixIt++;
  	    }
      }
    List<T>& operator<< (const T &t) 
      { P::push_back (t); 
      	return *this;
      }    
    List<T>& operator<< (T &&t) 
      { P::push_back (std::move (t)); 
      	return *this;
      }    
    template <typename U/*:<T>*/>
      List<T>& operator<< (const list<U> &other)
        { P::insert (P::end (), other. begin (), other. end ());
        	return *this;
        }
    template <typename U/*:<T>*/>
      List<T>& operator<< (const vector<U> &other)
        { P::insert (P::end (), other. begin (), other. end ());
        	return *this;
        }
    T popFront ()
      { if (P::empty ())
          throw runtime_error ("popFront() empty list");
        const T t = P::front ();
        P::pop_front ();
        return t;
      }
    T popBack ()
      { if (P::empty ())
          throw runtime_error ("popBack() empty list");
        const T t = P::back ();
        P::pop_back ();
        return t;
      }
  };



List<string> str2list (const string &s,
                       char c = ' ');
  // Invokes: findSplit(s,c)

string list2str (const List<string> &strList,
                 const string &sep = " ");



/* Usage:
  for (Iter<T> iter (t); iter. next (); )
    ...
*/
template <typename T>
  struct Iter final : Nocopy
  {
  private:
    T& t;
    typename T::iterator itNext;
    typename T::iterator it;
  public:
      

    explicit Iter (T &t_arg)
      : t (t_arg)
      , itNext (t_arg. begin ())
      {}

      
    bool next () 
      { it = itNext;
        if (it == t. end ())
          return false;
        itNext = it;
        itNext++;
        return true;
      }
    size_t getIndex () const
      { return (size_t) (it - t. begin ()); }
    typename T::value_type& operator* () const
      { return const_cast <typename T::value_type&> (*it); }  // *set::iterator = const set::value_type
    typename T::value_type* operator-> () const
      { return & const_cast <typename T::value_type&> (*it); }  // *set::iterator = const set::value_type
    typename T::value_type erase ()
      { typename T::value_type val = std::move (*it);
        itNext = t. erase (it); 
        return val;
      }
    void insert (const typename T::value_type &val)
      { itNext = t. insert (itNext, val) + 1; }
      // Opposite to erase()
  };



// File

constexpr char fileSlash = 
  #ifdef _MSC_VER
    '\\'
  #else  // UNIX
    '/'
  #endif
  ;

inline string getFileName (const string &path)  
  { const size_t pos = path. rfind ('/');
  	if (pos == string::npos)
  		return path;
  	return path. substr (pos + 1);
  }

inline string getDirName (const string &path)  
  { const size_t pos = path. rfind ('/');
  	if (pos == string::npos)
  		return noString;
  	return path. substr (0, pos + 1);
  }

inline bool isDirName (const string &path)
  { return isRight (path, "/"); }


inline void addDirSlash (string &dirName)
  { if (! dirName. empty () && ! isDirName (dirName))
    	dirName += "/";
  }

inline string shellQuote (string s)
  { replaceStr (s, "\'", "\'\"\'\"\'");
  	return "\'" + s + "\'";
  }

inline string getFileExtension (const string &path)  
  { const size_t pos = path. rfind ('.');
  	if (pos == string::npos)
  		return noString;
  	return path. substr (pos + 1);
  }

inline string trimExtension (const string &path)
  {
    const size_t pos = path. rfind ('.');
    if (pos == string::npos)
      return path;
    return path. substr (0, pos);
  }

bool fileExists (const string &fName);

inline void checkFile (const string &fName)
  { if (! fileExists (fName))
      throw runtime_error ("File " + strQuote (fName) + " does not exist");
  }

streamsize getFileSize (const string &fName);

void copyText (const string &inFName,
               size_t skipLines,
               ostream &os);


#ifndef _MSC_VER
  inline void moveFile (const string &from,
                        const string &to)
    { if (const int code = ::rename (from. c_str (), to. c_str ()))
        throw runtime_error ("Cannot move file + " + shellQuote (from) + " to " + shellQuote (to) + " (" + to_string (code) + ")");
    }

  inline void removeFile (const string &fName)
    { if (const int code = ::remove (fName. c_str ()))
        throw runtime_error ("Cannot remove file + " + shellQuote (fName) + " (" + to_string (code) + ")");
    }

      
  inline string path2canonical (const string &path)
    { if (char* p = realpath (path. c_str (), nullptr))
      { const string s (p);
        free (p);
        return s;
      }
      throw runtime_error ("path2canonical " + shellQuote (path));
    }
  
  
  enum class Filetype {none, dir, terminal, disk, file, pipe, link, socket};

  string filetype2name (Filetype t);

  Filetype getFiletype (const string &path,
                        bool expandLink);
#endif



// Directory

#ifndef _MSC_VER
  inline bool directoryExists (const string &dirName)
    { return getFiletype (dirName, true) == Filetype::dir; }
  
  void createDirectory (const string &dirName);

  inline void removeDirectory (const string &dirName)
    { std::filesystem::remove_all (dirName); }
    // With its contents
                        
  string makeTempDir ();
    // And test

  void concatTextDir (const string &inDirName,
                      const string &outFName);
#endif



struct Dir
{
  List<string> items;
    // Simplified: contains no redundant "", ".", ".."
    // items.front().empty (): root


  explicit Dir (const string &dirName);
  Dir () = default;

    
  string get () const
    { return items. empty () 
               ? string (1, '.') 
               : nvl (list2str (items, string (1, fileSlash)), string (1, fileSlash)); 
    }
  string getParent () const
    { if (items. empty ())
        return "..";
      if (items. size () == 1 && items. front (). empty ())
        throw runtime_error ("Cannot get the parent directory of the root");
      const Dir parent (get () + "/..");
      return parent. get ();
    }
#ifndef _MSC_VER
  size_t create ();
    // Return: number of directories created
    //         0 <=> complete directory exists
#endif
};



void setSymlink (string fromFName,
                 const string &toFName,
                 bool fromPathIsAbsolute);




/////////////////////////////////////// streams //////////////////////////////////////////

// Binary streams
// File content is platform-dependent

template <typename T>
  inline void writeBin (ostream &f,
                        const T &t)
    { f. write (reinterpret_cast <const char*> (& t), sizeof (t)); }

template <typename T>
  inline void readBin (istream &f,
                       T &t)
    { f. read (reinterpret_cast <char*> (& t), sizeof (t)); }

  

// istream

bool getChar (istream &is,
              char &c);
  // Output: c if (bool)Return
              
void skipLine (istream &is);

void readLine (istream &is,
               string &s);
  // Output: s



struct Istringstream final : istringstream
{
  Istringstream () = default;
  void reset (const string &s)
    { clear ();
      str (s);
    }
};



// ostream

bool isRedirected (const ostream &os);



struct Offset 
// Linked to the ostream os
// Not thread-safe
{
private:
	static size_t size;
public:
	static constexpr size_t delta {2};


	Offset ()
  	{ size += delta; }
 ~Offset () 
  	{ size -= delta; }


  static void newLn (ostream &os) 
    { os << endl << string (size, ' '); }
};



struct Color
{
  enum Type { none    =  0
              // UNIX
            , black   = 30  
            , red     = 31 	
            , green   = 32
            , yellow  = 33
            , blue    = 34 	
            , magenta = 35 	
            , cyan    = 36 	
            , white   = 37
            };
            
  static string code (Type color = none, 
                      bool bright = false)
    { 
		#ifdef _MSC_VER
		  return noString
		#else
    	return string ("\033[") + (bright ? "1;" : "") + to_string ((int) color) + "m"; 
		#endif
    }
};
  
  
inline string colorize_raw (const string &s,
                            Color::Type color,
                            bool bright,
                            bool screen)
  { if (! screen)
      return s;
    return Color::code (color, bright) + s + Color::code ();
  }

inline string colorize (const string &s,
                        bool screen)
  { return colorize_raw (s, Color::white, true, screen); }  
    
inline string colorizeUrl (const string &s,
                           bool screen)
  { return colorize_raw (s, Color::blue, true, screen); }  
  
inline string colorizeDir (const string &s,
                           bool screen)
  { return colorize_raw (s, Color::green, false, screen); }  


class OColor
// Output color
{
	ostream &o;
	const bool active;
public:
  OColor (ostream &o_arg,
          Color::Type color, 
          bool bright,
          bool active_arg)
	  : o (o_arg)
	  , active (active_arg && ! isRedirected (o_arg))
    { if (active) 
        o << Color::code (color, bright);
    }
 ~OColor ()    
    { if (active) 
        o << Color::code ();
    }
};



class ONumber
// Output number format
{
	ostream &o;
  const streamsize prec_old;
	const ios_base::fmtflags flags_old;
public:
	ONumber (ostream &o_arg,
	         streamsize precision,
	         bool scientific_arg)
	  : o (o_arg)
	  , prec_old (o. precision ())
	  , flags_old (o. flags ())
	  { if (scientific_arg)
	  	  o << scientific;
	  	else
	  		o << fixed;
      o. precision (precision);
	  }
 ~ONumber () 
    { o. flags (flags_old); 
      o. precision (prec_old); 
    }
};



struct COutErr : Singleton<COutErr>
{
  const bool both;
  
  COutErr ()
    : both (
           #ifdef _MSC_VER
             true
           #else
             sameFiles (fileno (stdout), fileno (stderr))
           #endif
           )
    {}
#ifndef _MSC_VER
private:
  static bool sameFiles (int fd1, 
                         int fd2);
public:
#endif
  template <typename T>
    const COutErr& operator<< (const T &val) const
      { cout << val;
        if (! both)
          cerr << val;
        return *this;
      }
  const COutErr& operator<< (ostream& (*pfun) (ostream&)) const
    { pfun (cout);
      if (! both)
        pfun (cerr);
      return *this;
    }
};



extern const COutErr couterr;



inline void section (const string &title,
                     bool useAlsoCout)
  { const Color::Type color = Color::cyan;
    const bool bright = false;
    if (useAlsoCout)
    {
      { const OColor oc1 (cout, color, bright, true);
        const OColor oc2 (cerr, color, bright, true);
        couterr << title /*<< " ..."*/;
      }
      couterr << endl;
    }
    else
    {
      { const OColor oc (cerr, color, bright, true);
        cerr << title /*<< " ..."*/;
      }
      cerr << endl;
    }
  }



template <typename T>
  inline void report (ostream &os,
			                const string &name,
			                const T &value)
    {	os << name << '\t' << value << endl; }
  


struct IFStream : ifstream
// Text file
{
  IFStream () = default;
	explicit IFStream (const string &pathName);
};



struct OFStream : ofstream
// Text file
{
	OFStream () = default;
	OFStream (const string &dirName,
	          const string &fileName,
	          const string &extension)
	  { open (dirName, fileName, extension); }
	explicit OFStream (const string &pathName)
	  { open (noString, pathName, noString); }
	  
	  
	static void create (const string &pathName)
	  { try { OFStream f (pathName); }
		    catch (...) { throw runtime_error ("Cannot create output file " + strQuote (pathName)); }
	  }
	static void prepare (const string &pathName)
	  { if (! pathName. empty ())
	    { create (pathName);
        removeFile (pathName);
      }
	  }
	void open (const string &dirName,
	           const string &fileName,
	           const string &extension);
	  // Input: !fileName.empty()
};



struct Cout : Singleton<Cout>
{
  unique_ptr<OFStream> f;
  
  
  explicit Cout (const string &fName)
    { if (! fName. empty ())
        f. reset (new OFStream (fName));
    }
    
  
  ostream& operator* ()
    { if (f)
        return *f;
      return cout;
    }
};



inline streamsize double2decimals (double r)
  { return r ? (streamsize) max<long> (0, (long) (ceil (- log10 (fabs (r)) + 1))) : 0; }




//////////////////////////////////////////////////////////////////////////////////////

void exec (const string &cmd,
           const string &logFName = noString);
  // Input: logFName: log file populated by cmd, to include into exception::what() if cmd fails

#ifndef _MSC_VER
  string which (const string &progName);
    // Return: isRight(,"/") or empty() if there is no path
#endif



struct Lock
{
  const bool active;
private:
  std::mutex& mtx;
public:
    
    
  explicit Lock (std::mutex &mtx_arg,
                 bool active_arg = true)
    : active (active_arg)
    , mtx (mtx_arg)
    { if (active)
        mtx_arg. lock (); 
    }
 ~Lock ()
    { if (active)
        mtx. unlock (); 
    }
};



// Threads

extern size_t threads_max;
  // >= 1
  
extern thread::id main_thread_id;
	
bool isMainThread ();

#ifndef _MSC_VER
  size_t get_threads_max_max ();
#endif



struct Threads : Singleton<Threads>
// Usage: { Threads th; th << ...; main_thread_process(); }
{
private:
	static size_t threadsToStart;
	vector<thread> threads;
	static bool quiet;
public:


  explicit Threads (size_t threadsToStart_arg,
                    bool quiet_arg = false);
 ~Threads ();

  	
	static bool empty () 
	  { return ! threadsToStart; }
	static bool isQuiet () 
	  { return quiet; }
	size_t getAvailable () const
	  { return threadsToStart < threads. size () ? 0 : (threadsToStart - threads. size ()); }
	Threads& operator<< (thread &&t)
	  { if (! getAvailable ())
	  	  throwf ("Too many threads created");
	  	try { threads. push_back (std::move (t)); }
	  	  catch (const exception &e) 
	  	    { throw runtime_error (string ("Cannot start thread\n") + e. what ()); }
	  	return *this;
	  }
	bool exec (const string cmd,
	           size_t cmdThreads = 1)
	  { if (cmdThreads && getAvailable () >= cmdThreads)
	    { *this << thread (Common_sp::exec, cmd, noString);
	      threadsToStart -= cmdThreads - 1;
	      return true;
	    }
      Common_sp::exec (cmd);
	    return false;
	  }
};



template <typename Func, typename Res, typename... Args>
  void arrayThreads (bool quiet,
                     const Func& func,
                     size_t i_max,
                     vector<Res> &results,
                     Args&&... args)
  // Input: void func (size_t from, size_t to, Res& res, Args...)
  // Optimial thread_num minimizes (Time_SingleCPU/thread_num + Time_openCloseThread * (thread_num - 1)), which is sqrt(Time_SingleCPU/Time_openCloseThread)
  {
  	if (threads_max < 1)
  	  throwf ("threads_max < 1");
		results. clear ();
		results. reserve (threads_max);
  	if (threads_max == 1 || i_max <= 1 || ! Threads::empty ())
  	{
  		results. push_back (Res ());
    	func (0, i_max, results. front (), std::forward<Args>(args)...);
  		return;
  	}
		size_t chunk = max<size_t> (1, i_max / threads_max);
		if (chunk * threads_max < i_max)
			chunk++;
		if (chunk * threads_max < i_max)
		  throwf ("chunk * threads_max < i_max");
		Threads th (threads_max - 1, quiet);
		for (size_t tn = 0; tn < threads_max; tn++)
	  {
	    const size_t from = tn * chunk;
	  	if (from >= i_max)
	  		break;
	    const size_t to = from + chunk;
	    results. push_back (Res ());
	    Res& res = results. back ();
	    if (to >= i_max)
	    {
	    	func (from, i_max, res, std::forward<Args>(args)...);
	    	break;
	    }
		  th << thread (func, from, to, ref (res), std::forward<Args>(args)...);
		}
  }




///////////////////////////////////////////////////////////////////////////////

template <typename T>
  struct Enumerate
  {
    vector<T> num2elem;
  private:
    unordered_map<T,size_t> elem2num;
  public:
    

    explicit Enumerate (size_t n)
      { num2elem. reserve (n);
        elem2num. rehash (n);
      }

      
    size_t size () const
      { return num2elem. size (); }
    size_t find (const T &t) const
      { if (const size_t* num = Common_sp::findPtr (elem2num, t))
          return *num;
        return no_index;
      }
    size_t add (const T &t)
      { size_t num = find (t);
        if (num == no_index)
        { num2elem. push_back (t);
          num = num2elem. size () - 1;
          elem2num [t] = num;
        }
        return num;
      }
  };



typedef  Enumerate<string>  Names;




//////////////////////////////////////////////////////////////////////////////////////

struct Xml
{
  struct File;


  struct Tag
  {
  private:
    File &f;
  public:
    const string name;
    const bool active;

    Tag (File &f_arg,
         const string &name_arg,
         bool active_arg = true);
   ~Tag ();
  };
  
  
  struct File 
  {
    friend Tag;
  protected:
  	bool isText {false};
  public:
    unique_ptr<const Tag> rootTag;

  protected:
  	File () = default;
  public:
    virtual ~File () 
      {}

  private:
    virtual void printRaw (const string &s) = 0;
    virtual void printRawText (const string &s) = 0;
    virtual void tagStart (const string &tag) = 0;
    virtual void tagEnd   (const string &tag) = 0;
  public:
  	void print (const string &s)
  	  { printRawText (s);
  	  	isText = true;
  	  }
		void text (const string &s) 
		  { const Tag xml (*this, "text");
		    print (s);
		 	}
    template <typename T>
      File& operator<< (const T &t)
        { ostringstream oss;
        	oss << t;
        	print (oss. str ());
          return *this;
        }
  };
  

  
  struct TextFile final : File
  // Tag::name: idenifier with possible '-'
  {
  private:
    struct XmlStream : OFStream
    { 
      explicit XmlStream (const string &pathName)
        : OFStream (pathName)
        { *this << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl; }  
    };
    XmlStream os;
  public:

    TextFile (const string &pathName,
		          const string &rootTagName)
      : os (pathName)
      { rootTag. reset (new Tag (*this, rootTagName)); }
   ~TextFile ()
      { rootTag. reset (); }

  private:
    void printRaw (const string &s) final
      { os << s; }
    void printRawText (const string &s) final;
    void tagStart (const string &tag) final;
    void tagEnd (const string &tag) final;
  };
  
  
  
  struct BinFile final : File
  // Binary XML
  //   <Data> ::= <nameIndex> <Data>* 0 0 <text> 0
  //     <nameIndex> ::= <byte> <byte>
  //   Number of different Tag::name's <= 2^16
  //   Tag::name: no '\0', '\n'
  {
  private:
  	ofstream os;
    Names names;
  public:

    BinFile (const string &pathName,
             const string &rootTagName)
      : os (pathName, ios_base::binary | ios_base::out)
      , names (10000)  // PAR
      { rootTag. reset (new Tag (*this, rootTagName)); }
   ~BinFile ();

  private:
    void printRaw (const string &s) final;
    void printRawText (const string &s) final
      { printRaw (s); }
    void tagStart (const string &tag) final;
    void tagEnd (const string &tag) final;
  };
};



extern unique_ptr<Xml::File> cxml;
  



//////////////////////////////////////// Root ///////////////////////////////////////////

struct Json;
struct JsonNull;
struct JsonInt;
struct JsonDouble;
struct JsonString;
struct JsonBoolean;
struct JsonContainer;
struct JsonArray;
struct JsonMap;  



struct Root
{
protected:
  Root () = default;
public:
  virtual ~Root () 
    {}
    // A desrtructor should be virtual to be automatically invoked by a descendant class destructor
  virtual Root* copy () const
    { throwf ("Root::copy() is not implemented"); }
    // Return: the same type    
  virtual void qc () const
    {}
    // Input: qc_on
  virtual void saveText (ostream& /*os*/) const 
    { throwf ("Root::saveText() is not implemented"); }
  void saveFile (const string &fName) const;
    // if fName.empty() then do nothing
    // Invokes: saveText()
  string str () const
    { ostringstream oss;
      saveText (oss);
      return oss. str ();
    }
  virtual void saveXml (Xml::File& /*f*/) const 
    { throwf ("Root::saveXml() is not implemented"); }
  virtual Json* toJson (JsonContainer* /*parent_arg*/,
                        const string& /*name_arg*/) const
    { throwf ("Root::toJson() is not implemented"); }
	virtual bool empty () const
	  { throwf ("Root::empty() is not implemented"); }
  virtual void clear ()
    { throwf ("Root::clear() is not implemented"); }
    // Postcondition: empty()
};

 

inline ostream& operator<< (ostream &os,
                            const Root &r) 
  { r. saveText (os);
    return os;
  }  



template <typename T /*Root*/> 
  struct AutoPtr final : unique_ptr<T>  
  {
  private:
  	typedef  unique_ptr<T>  P;
  public:
    
  
  	explicit AutoPtr (T* t = nullptr) 
  	  : P (t)
  	  {}
  	AutoPtr (const AutoPtr<T> &t) 
  	  : P (t. copy ())
  	  {}
  	AutoPtr (AutoPtr<T> &&t) 
  	  : P (t. release ())
  	  {}
  	AutoPtr<T>& operator= (T* t)
  	  { P::reset (t);
  	  	return *this;
  	  }
  	AutoPtr<T>& operator= (const AutoPtr<T> &t)
  	  { P::reset (t. copy ());
  	  	return *this;
  	  }
  	AutoPtr<T>& operator= (AutoPtr<T> &&t)
  	  { P::reset (t. release ());
  	  	return *this;
  	  }
  	T* copy () const
  	  { return P::get () ? static_cast<T*> (P::get () -> copy ()) : nullptr; }
  };



struct VirtNamed : Root
{
	virtual string getName () const = 0;
};


template <typename ParentValue/*:VirtNamed*/, typename Value/*:ParentValue*/>
  inline bool containsNamed (const map <string, const ParentValue*> &m,
                             const Value* t)
    { if (! t)
        return false;
      const auto& it = m. find (t->getName ());
      if (it == m. end ())
        return false;
      return it->second == t;
    }



struct Named : VirtNamed
{
  string name;
    // !empty(), no spaces at the ends, printable ASCII characeters


  Named () = default;
  explicit Named (const string &name_arg)
    : name (name_arg) 
    {}
  explicit Named (string &&name_arg)
    : name (std::move (name_arg))
    {}
  Named* copy () const override
    { return new Named (*this); } 
	string getName () const override
	  { return name; }


  void qc () const override;
  void saveText (ostream& os) const override
    { os << name; }
	bool empty () const override
	  { return name. empty (); }
  void clear () override
    { name. clear (); }
	static bool lessPtr (const Named* x,
	                     const Named* y)
	  { return x->name < y->name; }
};


inline string named2name (const Named* n)
  { return n ? n->name : "<anonymous>"; }




///////////////////////////// Container extensions ////////////////////////////////////////////

template <typename T>
  struct Vector : vector<T>
  {
  private:
  	typedef  vector<T>  P;
  public:
    ebool ascending {enull};
  	

  	Vector () = default;
  	explicit Vector (size_t n, 
  	                 const T &value = T ())
  	  : P (n, value)
  	  {}
    explicit Vector (initializer_list<T> init)
      : P (init)
      {}
  	explicit Vector (const vector<T> &other) 
      : P (other)
      {}
  	
  	
  private:
  	void checkIndex (const string &operation,
  	                 size_t index) const
  	  {
  	  #ifndef NDEBUG
  	    if (index >= P::size ())
  	      throw runtime_error ("Vector " + operation + " operation: index = " + to_string (index) + ", but size = " + to_string (P::size ()));
  	  #endif
  	  }
  public:
    void unsetSearchSorted ()
      { if (P::size () > 1)
          ascending = enull;
      }
  	typename P::reference operator[] (size_t index)
  	  { checkIndex ("assignment", index);
  	    return P::operator[] (index);
  	  }
  	typename P::const_reference operator[] (size_t index) const
  	  { checkIndex ("get", index);
  	    return P::operator[] (index);
  	  }
    void reserveInc (size_t inc)
      { P::reserve (P::size () + inc); }
    bool find (const T &value,
               size_t &index) const
  	  // Output: index: valid if (bool)Return
  	  { for (index = 0; index < P::size (); index++)
  	      if ((*this) [index] == value)
  	        return true;
  	    return false;
  	  }
    typename Vector<T>::const_iterator constFind (const T &value) const
  	  { typename Vector<T>::const_iterator it = P::begin (); 
  	    while (it != P::end ()) 
  	      if (*it == value)
  	        break;
  	      else
  	        it++;
  	    return it;
  	  }
    typename Vector<T>::iterator find (const T &value)
  	  { typename Vector<T>::iterator it = P::begin (); 
  	    while (it != P::end ()) 
  	      if (*it == value)
  	        break;
  	      else
  	        it++;
  	    return it;
  	  }
    bool contains (const T &value) const
      { return constFind (value) != P::end (); }
    bool containsAll (const vector<T> &other) const
      { for (const T& t : other)
          if (! contains (t))
            return false;
        return true; 
      }
    size_t indexOf (const T &value) const
      { size_t n = 0;
        for (const T& t : *this)
          if (value == t)
            return n;
          else
            n++;
        return no_index;
      }
    size_t countValue (const T &value) const
      { size_t n = 0;
        for (const T& t : *this)
          if (value == t)
            n++;
        return n;
      }
    bool overlapStart (size_t start,
                       const vector<T> &other) const
      { const size_t end = min (P::size (), start + other. size ());
        for (size_t i = start, j = 0; i < end; i++, j++)
          if ((*this) [i] != other [j])
            return false;
        return true;
      }
    size_t overlapStart_min (const vector<T> &other,
                             size_t start_min) const
      { for (size_t start = start_min; start < P::size (); start++)
          if (overlapStart (start, other))
            return start;
        return P::size ();
      }
      // a = v.overlapStart_min(v,1), a < v.size(), v.size() % a == 0 => a is the period of v
    void checkSorted () const
      { if (ascending == enull)
      	  throw runtime_error ("Vector is not sorted for search");
      }
    template <typename U>
      void checkSorted (const Vector<U> &other) const
        { checkSorted ();
          other. checkSorted ();
          if (ascending != other. ascending)
            throw logic_error ("Different ordering direction of two vectors");
        }
    Vector<T>& operator<< (const T &value)
      { P::push_back (value);
        unsetSearchSorted ();
      	return *this;
      }
    Vector<T>& operator<< (T &&value)
      { P::push_back (std::move (value));
        unsetSearchSorted ();
      	return *this;
      }
    template <typename U/*:<T>*/>
      Vector<T>& operator<< (const vector<U> &other)
        { reserveInc (other. size ());
          P::insert (P::end (), other. begin (), other. end ());
          unsetSearchSorted ();
        	return *this;
        }
    template <typename U/*:<T>*/>
      Vector<T>& operator<< (vector<U> &&other)
        { reserveInc (other. size ());
          for (U& t : other)
            P::push_back (std::move (t));
          unsetSearchSorted ();
          other. clear ();
        	return *this;
        }
    template <typename U/*:<T>*/>
      Vector<T>& operator<< (const list<U> &other)
        { reserveInc (other. size ());
          P::insert (P::end (), other. begin (), other. end ());
          unsetSearchSorted ();
        	return *this;
        }
    void setAll (const T &value)
      { for (T &t : *this)
      	  t = value;
      }
    void setAt (size_t index,
                T value)
      { if (index >= P::size ())
          P::resize (index + 1);
        if ((*this) [index] == T ())
          (*this) [index] = value;
        else
          throw runtime_error ("vector [" + to_string (index) +"] is not empty");
      }
    void eraseAt (size_t index)
      { eraseMany (index, index + 1); }
    void eraseMany (size_t from,
                    size_t to)
      { if (to < from)
          throw runtime_error ("Vector::eraseMany(): to < from");
        if (to == from)
          return;
        checkIndex ("eraseMany", to - 1);
        auto it1 = P::begin ();
        std::advance (it1, from);
        auto it2 = it1;
        std::advance (it2, to - from);
        P::erase (it1, it2);
      }
    void wipe ()
      { P::clear ();
      	P::shrink_to_fit ();
      }
    void reverse ()
      { for (size_t i = 0; i < P::size (); i++)
      	{ const size_t j = P::size () - 1 - i;
      		if (i >= j)
      			break;
      	  std::swap ((*this) [i], (*this) [j]);
      	}
        toggle (ascending);
      }
    const T& getRandom (Rand &rand) const
      { return (*this) [rand. get (P::size ())]; }
    void randomOrder ()
  		{ Rand rand (seed_global);
  			for (T &t : *this)
  	      std::swap (t, (*this) [(size_t) rand. get ((ulong) P::size ())]);
      	unsetSearchSorted ();
  		}
    void pop_back ()
      { 
      #ifndef NDEBUG
        if (P::empty ())
          throw runtime_error ("Empty vector pop_back");
      #endif
        P::pop_back ();
      }
    T pop (size_t n = 1)
      { T t = T ();
        while (n)
        { if (P::empty ())
            throw runtime_error ("Cannot pop an empty vector");
          t = (*this) [P::size () - 1];
      	  P::pop_back ();
          n--;
        }
      	return t;
      }
    template <typename Condition /*on T*/>
      size_t count (Condition cond) const
        { size_t n = 0;
          for (const T& t : *this)
            if (cond (t))
              n++;
          return n;
        }
    template <typename Condition /*on index*/>
      void filterIndex (const Condition cond)
        { size_t toDelete = 0;
          for (size_t i = 0, end_ = P::size (); i < end_; i++)
          { const size_t j = i - toDelete;
            if (j != i)
              (*this) [j] = std::move ((*this) [i]);
            if (cond (j))
              toDelete++;
          }
          while (toDelete)
          { P::pop_back ();
            toDelete--;
          }
        }
    template <typename Condition /*on value*/>
      void filterValue (const Condition cond)
        { size_t toDelete = 0;
          for (size_t i = 0, end_ = P::size (); i < end_; i++)
          { const size_t j = i - toDelete;
            if (j != i)
              (*this) [j] = std::move ((*this) [i]);
            if (cond ((*this) [j]))
              toDelete++;
          }
          while (toDelete)
          { P::pop_back ();
            toDelete--;
          }
        }

    void sort ()
      { if (ascending == etrue)
          return;
        Common_sp::sort (*this); 
        ascending = etrue;
      }
    template <typename StrictlyLess>
      void sort (const StrictlyLess &strictlyLess)
        { Common_sp::sort (*this, strictlyLess); 
          unsetSearchSorted ();
        }    
    bool lt (const T& t1,
             const T& t2) const
      { if (ascending == etrue)
          return t1 < t2;
        return t2 < t1;
      }
      // Requires: acsending != enull
    void sortBubble ()  
      // Fast if *this is almost sort()'ed
      { if (ascending == etrue)
          return;
        for (size_t i = 1; i < P::size (); i++)
  		    for (size_t j = i; j-- > 0;)
  		      if (lt ((*this) [j + 1], (*this) [j]))
          	  std::swap ((*this) [j + 1], (*this) [j]);
  		      else
  		      	break;
      	ascending = etrue;
      }

    size_t binSearch (const T &value,
                      bool exact = true) const
      // Return: if exact then no_index or vec[Return] = value else min {i : vec[i] >= value}
      { if (P::empty ())
      	  return no_index;
      	checkSorted ();
      	if (P::size () < 8)  // PAR
      	{
      	  for (size_t i = 0; i < P::size (); i++)
      	  {
      	    if (value == (*this) [i])
          	  return i;
      	    if (lt (value, (*this) [i]))
          	  return exact ? no_index : i;
          }
          return no_index;
      	}
      	size_t lo = 0;  // vec.at(lo) <= value
      	size_t hi = P::size () - 1;  
      	// lo <= hi
      	if (lt (value, (*this) [lo]))
      	  return exact ? no_index : lo;
      	if (lt ((*this) [hi], value))
      	  return no_index;
      	// at(lo) <= value <= at(hi)
      	for (;;)
      	{
  	    	const size_t m = (lo + hi) / 2;
  	    	if (       (*this) [m] == value
  	    		  || lt ((*this) [m],   value)
  	    		 )
  	    		if (lo == m)  // hi in {lo, lo + 1}
  	    			break;
  	    		else
  	    		  lo = m;
  	      else
  	      	hi = m;
  	    }
  	    if ((*this) [lo] == value)
  	    	return lo;
  	    if (! exact || (*this) [hi] == value)
  	    	return hi;
  	    return no_index;
      }
    template <typename U /* : T */>
      bool containsFast (const U &value) const
        { return binSearch (value) != no_index; }
    template <typename U /* : T */>
      bool containsFastAll (const Vector<U> &other) const
        { if (other. size () > P::size ())
      	    return false;
          for (const U& u : other)
            if (! containsFast (u))
              return false;
          return true;
        }
    template <typename U /* : T */>
      bool containsFastAll (const set<U> &other) const
        { if (other. size () > P::size ())
      	    return false;
          for (const U& u : other)
            if (! containsFast (u))
              return false;
          return true;
        }
    template <typename U /* : T */>
      bool containsFastAll (const unordered_set<U> &other) const
        { if (other. size () > P::size ())
      	    return false;
          for (const U& u : other)
            if (! containsFast (u))
              return false;
          return true;
        }
    template <typename U /* : T */, typename V>
      bool containsFastAll (const unordered_map<U,V> &other) const
        { if (other. size () > P::size ())
      	    return false;
          for (const auto& it : other)
            if (! containsFast (it. first))
              return false;
          return true;
        }
    template <typename U /* : T */>
      bool intersectsFast (const Vector<U> &other) const
        { if (P::empty ())
            return false;
          for (const U& u : other)
            if (containsFast (u))
              return true;
          return false;
        }
    template <typename U>
      bool intersectsFast_merge (const Vector<U> &other) const
        { checkSorted (other);
        	size_t i = 0;
        	const size_t otherSize = other. size ();
          for (const T& t : *this)
          { while (i < otherSize && lt (other [i], t))
              i++;
            if (i == otherSize)
              return false;
            if (other [i] == t)
              return true;
          }
          return false;
        }
    template <typename U /* : T */>
      bool intersects (const set<U> &other) const
        { if (other. empty ())
            return false;
          for (const T& t : *this)
            if (other. find (t) != other. end ())
              return true;
          return false;
        }
    template <typename U /* : T */>
      void setMinus (const Vector<U> &other)
        { filterIndex ([&] (size_t i) { return other. containsFast ((*this) [i]); }); }
        
    size_t findDuplicate (bool checkSortedP = true) const
      { if (P::size () <= 1)
          return no_index;
        if (checkSortedP)
          checkSorted ();
        for (size_t i = 1; i < P::size (); i++)
          if ((*this) [i] == (*this) [i - 1])
            return i;
        return no_index;
      }
    bool isUniq (bool checkSortedP = true) const
      { return findDuplicate (checkSortedP) == no_index; }
    template <typename Equal /*bool equal (const T &a, const T &b)*/>
  	  void uniq (const Equal &equal)
  	    { if (P::size () <= 1)
  	        return;
  	      size_t j = 1;  
  	      for (size_t i = 1; i < P::size (); i++)
  	        if (! equal ((*this) [i], (*this) [i - 1]))
  	        { if (j != i)
  	            (*this) [j] = (*this) [i];
  	          j++;
  	        }
  	      P::resize (j);
  	    }
  	void uniq ()
  	  { uniq ([] (const T& a, const T& b) { return a == b; }); }
  	bool checkUniq () const
  	  { Vector<T> other (*this);
  	    other. sort ();
  	    return other. isUniq ();
  	  }
    size_t getIntersectionSize (const Vector<T> &other) const
      // Input: *this, vec: unique
      { if (other. empty ())
          return 0;
        checkSorted (other);
        size_t n = 0;
        size_t j = 0;
        for (const T& x : *this)
        { while (lt (other [j], x))
          { j++;
            if (j == other. size ())
              return n;
          }
          if (other [j] == x)
            n++;
        }
        return n;
      }
    Vector<T> getIntersection (const Vector<T> &other) const
      // Input: *this, vec: unique
      { Vector<T> res;
        if (other. empty ())
          return res;
        checkSorted (other);
        size_t j = 0;
        for (const T& x : *this)
        { while (lt (other [j], x))
          { j++;
            if (j == other. size ())
              return res;
          }
          if (other [j] == x)
            res << x;
        }
        res. ascending = ascending;
        return res;
      }
    bool setIntersection (const Vector<T> &other)
      { other. checkSorted ();      
        if (P::empty ())
          *this = other;
        else
        { Vector<T> vec (getIntersection (other));
          *this = std::move (vec);
        }
        ascending = other. ascending;
        return ! P::empty ();
      }
    Vector<T> getUnion (const Vector<T> &other) const
      // Input: *this, vec: unique
      { Vector<T> res;  res. reserve (P::size () + other. size ());      
        if (other. empty ())
          return *this;
        if (P::empty ())
          return other;
        checkSorted (other);
        size_t j = 0;
        for (const T& x : *this)
        { while (j < other. size () && lt (other [j], x))
          { res << other [j];
            j++;
          }
          if (j < other. size () && other [j] == x)
            j++;
          res << x;
        }
        while (j < other. size ())
        { res << other [j];
          j++;
        }
        res. ascending = ascending;
        return res;
      }
    void setUnion (const Vector<T> &other)
      { if (P::empty ())
          *this = other;
        else
        { Vector<T> vec (getUnion (other));
          *this = std::move (vec);
        }
      }

    bool operator< (const Vector<T> &other) const
      // Lexicographic comparison
      { const size_t n = std::min (P::size (), other. size ());
        for (size_t i = 0; i < n; i++)
      	{ if ((*this) [i] < other [i]) return true;
      		if (other [i] < (*this) [i]) return false;
        }
        return P::size () < other. size ();
      }
  };



// Vector<bool>

inline bool exists (const Vector<bool> &vec) 
  { for (const bool b : vec)
      if (b)
        return true;
    return false;
  }

inline size_t count (const Vector<bool> &vec) 
  { size_t n = 0;
    for (const bool b : vec)
      if (b)
        n++;
    return n;
  }



template <typename T>
  Vector<T> diff2vec (const unordered_set<T> &a,
                      const unordered_set<T> &b)
  { Vector<T> v;  v. reserve (a. size ());
    for (const T& t : a)
      if (! contains (b, t))
        v << t;
    return v;
  }  


template <typename Key /*VirtNamed*/, typename Value>
  Vector<pair<string,const Value*>> map2sortedVec (const map <const Key*, Value> &m)
    { Vector<pair<string,const Value*>> vec;  vec. reserve (m. size ());
	    for (const auto& it : m)
	    { assert (it. first);
	    	assert (& it. second);
	    	vec << pair<string,const Value*> (it. first->getName (), & it. second);
	    }
	    vec. sort ();    	
	    return vec;
    }
    


template <typename T /* : Root */>
  struct VectorPtr : Vector <const T*>
  {
  private:
  	typedef  Vector <const T*>  P;
  public:


    VectorPtr () = default;
  	explicit VectorPtr (size_t n, 
  	                    const T* value = nullptr)
  	  : P (n, value)
  	  {}
    explicit VectorPtr (initializer_list<const T*> init)
      : P (init)
      {}
  	template <typename U>
    	explicit VectorPtr (const vector<const U*> &other)
    	  : P ()
    	  { P::reserve (other. size ());
    	    insertAll (*this, other);  	    
    	  }	  
  	template <typename U>
    	explicit VectorPtr (const Vector<const U*> &other)
    	  : P ()
    	  { P::reserve (other. size ());
    	    insertAll (*this, other);
    	    P::ascending = other. ascending;
    	  }	  
  	template <typename U>
    	explicit VectorPtr (const list<const U*> &other)
    	  : P ()
    	  { P::reserve (other. size ());
    	    insertAll (*this, other);
    	  }	  
    template <typename U>
    	VectorPtr<T>& operator= (const Vector<const U*> &other)
    	  { P::clear ();
    	    P::reserve (other. size ());
    	    insertAll (*this, other);
    	    P::ascending = other. ascending;
    	    return *this;
    	  }	  


    VectorPtr<T>& operator<< (const T* value)
      { P::operator<< (value); 
      	return *this;
      }
    template <typename U/*:<T>*/>
      VectorPtr<T>& operator<< (const VectorPtr<U> &other)
        { P::operator<< (other); 
          return *this;
        }
    template <typename U/*:<T>*/>
      VectorPtr<T>& operator<< (VectorPtr<U> &&other)
        { P::operator<< (std::move (other)); 
          return *this;
        }
  	void deleteData ()
  	  {	for (const T* t : *this)
			    delete t;
  			P::clear ();
  	  }
    void erasePtr (size_t index)
      { delete (*this) [index];
        P::eraseAt (index);
      }
    void sortPtr ()
      { P::sort (lessPtr<T>); 
  		  P::unsetSearchSorted ();
      }
    void sortBubblePtr ()
      { for (size_t i = 1; i < P::size (); i++)
  		    for (size_t j = i; j-- > 0;)
  		      if (* (*this) [j + 1] < * (*this) [j])
          	  std::swap ((*this) [j], (*this) [j + 1]);
  		      else
  		      	break;
  		  P::unsetSearchSorted ();
      }
  };



template <typename T>
  bool equalPtr (const VectorPtr<T> &a,
                 const VectorPtr<T> &b)
    { if (a. size () != b. size ())
        return false;
      for (size_t i = 0; i < a. size (); i++)
        if (! (* a [i] == * b [i]))
          return false;
      return true;
    }

template <typename T>
  ebool compPtr (const VectorPtr<T> &a,
                 const VectorPtr<T> &b)
    { if (a. size () < b. size ())
        return etrue;
      if (a. size () > b. size ())
        return efalse;
      for (size_t i = 0; i < a. size (); i++)
      {
        if (* a [i] < * b [i])
          return etrue;
        if (* b [i] < * a [i])
          return efalse;
      }
      return enull;
    }


template <typename Key /*VirtNamed*/>
  Vector<pair<string,const Key*>> vec2sorted (const VectorPtr<Key> &in)
    { Vector<pair<string,const Key*>> vec;  vec. reserve (in. size ());
	    for (const Key* key : in)
	    { assert (key);
	    	vec << pair<string,const Key*> (key->getName (), key);
	    }
	    vec. sort ();    	
	    return vec;
    }



template <typename T /* : Root */>
  struct VectorOwn : VectorPtr<T>
  {
  private:
  	typedef  VectorPtr<T>  P;
  public:


    VectorOwn () = default;
  	VectorOwn (const VectorOwn<T> &other) 
  	  : P ()
  	  { *this = other; } 
  	VectorOwn<T>& operator= (const VectorOwn<T> &other) 
  	  { P::deleteData ();
  	  	P::reserve (other. size ());
  	  	for (const T* t : other)
  	  	  P::push_back (static_cast <const T*> (t->copy ()));
  	  	P::ascending = other. ascending;
  	  	return *this;
  	  }
  	VectorOwn (VectorOwn<T> &&other)
  	  : P ()
  	  { *this = std::move (other); }
  	VectorOwn<T>& operator= (VectorOwn<T> &&other)
  	  { P::deleteData ();
  	    P::operator= (std::move (other)); 
  	  	P::ascending = other. ascending;
  	    return *this;
  	  }
  	explicit VectorOwn (const VectorPtr<T> &other)
  	  : P ()
  	  { *this = other; }
  	VectorOwn<T>& operator= (const VectorPtr<T> &other) 
  	  { P::deleteData ();
  	    P::operator= (other); 
  	    P::ascending = other. ascending;
  	    return *this;
  	  }
  	explicit VectorOwn (VectorPtr<T> &&other)
  	  : P ()
  	  { *this = std::move (other); }
  	VectorOwn<T>& operator= (VectorPtr<T> &&other) 
  	  { P::deleteData ();
  	    P::operator= (std::move (other)); 
  	  	P::ascending = other. ascending;
  	    return *this;
  	  }
   ~VectorOwn ()
      { 
      #pragma GCC diagnostic ignored "-Wnull-dereference"
        P::deleteData (); 
      #pragma GCC diagnostic warning "-Wnull-dereference"
      }


    VectorOwn<T>& operator<< (const T* value)
      { P::operator<< (value); 
      	return *this;
      }
    template <typename U/*:<T>*/>
      VectorOwn<T>& operator<< (const VectorPtr<U> &other)
        { P::operator<< (other); 
          return *this;
        }


    class Stack
    {
      VectorOwn<T> &vec;
      const size_t size;
      
    public:
      explicit Stack (VectorOwn<T> &vec_arg)
        : vec (vec_arg)
        , size (vec_arg. size ())
        {}
     ~Stack ()
        { while (vec. size () > size)
          { const T* t = vec. pop ();
            delete t;
          }
        }
        
      void release (VectorPtr<T> &out)
        { out. reserve (vec. size () - size);
          while (vec. size () > size)
            out << vec. pop ();
        }
    };
  };



template <typename T>
  struct Set : set<T>
  {
  private:
  	typedef  set<T>  P;
  public:
  	bool universal {false};
  	  // true => empty()
  	

  	Set () = default;
  	explicit Set (bool universal_arg)
  	  : universal (universal_arg)
  	  {}
  	template <typename U, typename V>
  	  explicit Set (const map<U,V> &other)
  	    { operator= (other); }
  	template <typename U, typename V>
  	  explicit Set (const unordered_map<U,V> &other)
  	    { operator= (other); }
  	template <typename U>
  	  explicit Set (const vector<U> &other)
  	    { operator= (other); }
  	template <typename U, typename V>
    	Set<T>& operator= (const map<U,V> &other)
    	  { P::clear ();
    	    universal = false;
    	    for (const auto& it : other)
  	        P::insert (it. first);
  	      return *this;
    	  }
  	template <typename U, typename V>
    	Set<T>& operator= (const unordered_map<U,V> &other)
    	  { P::clear ();
    	    universal = false;
    	    for (const auto& it : other)
  	        P::insert (it. first);
  	      return *this;
    	  }
  	template <typename U>
    	Set<T>& operator= (const vector<U> &other)
    	  { P::clear ();
    	    universal = false;
    	    for (const U& u : other)
  	        P::insert (u);
  	      return *this;
    	  }
    bool operator== (const Set<T> &other) const
      { return universal
                 ? other. universal ? true : false
                 : other. universal
                   ? false
                   :    P::size () == other. size ()
                     && containsAll (other);
      }


    bool contains (const T& el) const
      { return universal || P::find (el) != P::end (); }
    T front () const
      { return * P::begin (); }
    T back () const
      { return * std::prev (P::end ()); }
    template <typename U /* : T */>
      bool containsAll (const Set<U> &other) const
        { if (universal)
      	    return true;
      	  if (other. universal)
      		  return false;
      	  if (other. size () > P::size ())
      	    return false;
          for (const U& u : other)
            if (! contains (u))
              return false;
          return true;
        }
    bool intersects (const Set<T> &other) const
       { if (universal && other. universal)
       	   return true;
       	 if (universal)
       	   return ! other. empty ();
       	 if (other. universal)
       	   return ! P::empty ();
       	 if (P::size () < other. size ())
           return intersects_ (other);
         return other. intersects_ (*this);
       }
  private:
    bool intersects_ (const Set<T> &other) const
      { if (other. empty ())
          return false;
        for (const T& t : *this)
          if (other. contains (t))
            return true;
        return false;
      }
  public:
    bool intersects (const unordered_set<T> &other) const
       { if (universal)
       	   return ! other. empty ();
         if (other. empty ())
           return false;
         for (const T& t : *this)
           if (contains (other, t))
             return true;
         return false;     
       }

    Set<T>& operator<< (const T &el)
      { if (! universal)
      	  P::insert (el);  // Slow
      	return *this;
      }
    Set<T>& operator<< (const Set<T> &other)
      { if (! universal)
      	{ if (other. universal)
      	  { P::clear ();
      	  	universal = true;
      	  }
      	  else
      	    P::insert (other. begin (), other. end ());
      	}
      	return *this;
      }
    template <typename From>
      void insertAll (const From &from)
        { P::insert (from. begin (), from. end ()); }
    bool addUnique (const T& el)
      { if (contains (el))
          return false;
        operator<< (el);
        return true;
      }
  	void intersect (const Set<T> &other) 
  		{ if (other. universal)
  			  return;
  			if (universal)
  			{ operator= (other);
  				return;
  			}
        for (Iter<Set<T>> iter (*this); iter. next (); )
  				if (! other. contains (*iter))
  					iter. erase ();
  		}
  	void intersect (const Vector<T> &other) 
  		{ if (universal)
  			{ operator= (other);
  				return;
  			}
        for (Iter<Set<T>> iter (*this); iter. next (); )
  				if (! other. containsFast (*iter))
  					iter. erase ();
  		}
  	template <typename U, typename V>
    	void intersect (const map<U,V> &other) 
    		{ if (universal)
    			{ operator= (other);
    				return;
    			}
          for (Iter<Set<T>> iter (*this); iter. next (); )
    				if (! Common_sp::contains (other, *iter))
    					iter. erase ();
    		}
  	template <typename U, typename V>
    	void intersect (const unordered_map<U,V> &other) 
    		{ if (universal)
    			{ operator= (other);
    				return;
    			}
          for (Iter<Set<T>> iter (*this); iter. next (); )
    				if (! Common_sp::contains (other, *iter))
    					iter. erase ();
    		}
  	size_t intersectSize (const Set<T> &other) const
  	  // Return: universal <=> SIZE_MAX
  		{ if (other. universal)
  			  return universal ? SIZE_MAX : P::size ();
  			if (universal)
  				return other. size ();
  		  size_t n = 0;
  		  if (! other. empty ())
    		  for (const T& t : *this)
    				if (other. contains (t))
    					n++;
  			return n;
  		}
    size_t setMinus (const Set<T> &other)
      { if (universal)
          throw runtime_error ("setMinus:universal");
      	size_t n = 0;
      	if (other. universal)
      	{ n = P::size ();
      		P::clear ();
      	}
      	else
  	    	for (const T& t : other)
  	        n += P::erase (t);
        return n;
      }
  };



template <typename T>
  void setMove (Set<T>* from,
	              Set<T>* to,
	              T el)
    { if (from == to)
    	  return;
    	assert (! from || ! from->universal);
    	assert (! to   || ! to  ->universal);
    	if (from)
    	  from->erase (el); 
    	if (to)
    	  to->insert (el);
    }



template <typename T, typename U /* : T */>
  bool containsFastAll (const Vector<T> &vec,
                        const Set<U> &other)
    { if (other. universal)
  		  return false;
  	  if (other. size () > vec. size ())
  	    return false;
      for (const U& u : other)
        if (! vec. containsFast (u))
          return false;
      return true;
    }



struct StringVector : Vector<string>
{
private:
	typedef  Vector<string>  P;
public:
	

  StringVector () = default;
  explicit StringVector (initializer_list<string> init)
    : P (init)
    {}
  explicit StringVector (const Set<string> &from)
    { insertAll (*this, from); 
      ascending = etrue;
    }
  StringVector (const string &fName,
                size_t reserve_size,
                bool trimP);
  StringVector (const string &s, 
                char sep,
                bool trimP);
  explicit StringVector (size_t n)
    : P (n, noString)
    {}


  string toString (const string& sep = noString) const;
  bool same (const StringVector &vec,
             const Vector<size_t> &indexes) const;
  void to_xml (Xml::File &f,
               const string &tag);
    // XML: <tag> <item>at(0)</item> <item>at(1)</item> ... </tag>
    // Invokes: sort(), clear()


  struct Hasher 
  {
    size_t operator () (const StringVector& vec) const 
    { size_t ret = 0;
      for (const string& s : vec) 
        ret ^= hash<string>() (s);
      return ret;
    }
  };
};



template <typename Key /*VirtNamed*/>
  StringVector set2vec (const set<const Key*> &s)
    { StringVector vec;  vec. reserve (s. size ());
	    for (const Key* key : s)
	    { assert (key);
	    	vec << key->getName ();
	    }
	    return vec;
    }


inline void aggregate (string& aggregation,
                       const string &addition,
                       char aggr_sep)
  { StringVector vec (aggregation, aggr_sep, true);
    vec << addition;
    vec. sort ();
    vec. uniq ();
    aggregation = vec. toString (string (1, aggr_sep));
  }



struct Csv 
// Line of Excel .csv-file
{
private:
  const string &s;
  size_t pos {0};
public:

  
  explicit Csv (const string &s_arg)
    : s (s_arg)
    {}
  
  
  bool goodPos () const
    { return pos < s. size (); }
  string getWord ();
    // Return: Next word
    // Requires: goodPos()
private:
  void findChar (char c)
    { while (goodPos () && s [pos] != c)
        pos++;
    }
};

  
  
StringVector csvLine2vec (const string &line);
  // Invokes: Csv



template <typename T>  
  struct RandomSet
  // Set stored in a vector for a random access
  {
  private:
    Vector<T> vec;
    unordered_map<T,size_t/*index in vec*/> um;
    // Same elements
  public:
    

    RandomSet () = default;
    void reset (size_t num)
      { vec. clear ();  vec. reserve (num);
        um.  clear ();  um. rehash (num);
      }
    void qc () const
      { if (! qc_on)
          return;
        if (vec. size () != um. size ())
          throwf ("RandomSet: qc");
      }

      
    // Time: O(1)
    bool empty () const
      { return vec. empty (); }
    size_t size () const
      { return vec. size (); }
    bool insert (const T &t)
      { if (contains (um, t))
          return false;
        um [t] = vec. size ();;
        vec << t;
        return true;
      }
    bool erase (const T &t)
      { if (vec. empty ())
          return false;
        const auto it = um. find (t);
        if (it == um. end ())
          return false;
        um [vec. back ()] = it->second;
        vec [it->second] = vec. back ();
        vec. pop_back ();
        um. erase (t);
        return true;
      }
    const Vector<T>& getVec () const
      { return vec; }
  };
  


template <typename T>
  struct Heap 
  // Priority queue
  // Heap property: comp(arr[parent(index)],arr[index]) >= 0
  // More operations than in std::priority_queue
  {
  private:
    Vector<T*> arr;
      // Elements are not owned by arr
    CompareInt comp {nullptr};
      // !nullptr
    typedef void (*SetHeapIndex) (T &item, size_t index);
      // Example: item.heapIndex = index
    SetHeapIndex setHeapIndex {nullptr};
      // Needed to invoke increaseKey()
  public:


    explicit Heap (const CompareInt &comp_arg,
  					       const SetHeapIndex &setHeapIndex_arg = nullptr,
  					       size_t toReserve = 0)
      : comp (comp_arg)
      , setHeapIndex (setHeapIndex_arg)
      { arr. reserve (toReserve); }
  private:
    static void throwError (const string &str) 
      { throw runtime_error ("Heap: " + str); }
  public:


    bool empty () const 
      { return arr. empty (); }
    size_t size () const
      { return arr. size (); }
    Heap& operator<< (T* item)
      { if (! item)
          throwError ("null item");
        arr << item;
        increaseKey (arr. size () - 1);
        return *this;
      }
    void increaseKey (size_t index)
      { T* item = arr [index];
        size_t p = no_index;
        while (index && comp (arr [p = parent (index)], item) < 0)
        { assign (arr [p], index);
          index = p;
        }
        assign (item, index);
      }
    void decreaseKey (size_t index)
      { heapify (index, arr. size ()); }
    T* getMaximum () const
      { if (arr. empty ()) 
      	  throwError ("getMaximum");
        return arr [0];
      }
    void deleteMaximum ()
      // Time: O(1) amortized
      { if (arr. empty ()) 
      	  throwError ("deleteMaximum");
        T* item = arr. back ();
        arr. pop_back ();
        if (arr. empty ())
          return;
        assign (item, 0);
        reinsertMaximum ();
      }
    void reinsertMaximum ()
      // Do this if the key of getMaximum() is changed
      { heapify (0, arr. size ()); }
    void sort ()
      { if (arr. empty ())
          return;
        for (size_t i = arr. size () - 1; i > 0; i--)
        { swap (0, i);
          heapify (0, i);
        }
      }
  private:
    static size_t parent (size_t index) 
      { if (! index)  throwError ("parent");
        return (index + 1) / 2 - 1;
      }
    static size_t left (size_t index) 
      { return 2 * (index + 1) - 1; }
    static size_t right (size_t index) 
      { return left (index) + 1; }
    void assign (T* item,
                 size_t index)
      { arr [index] = item;
        if (setHeapIndex)
          setHeapIndex (*item, index);
      }
    void swap (size_t i,
               size_t j)
      { if (i == j)
          return;
        T* item = arr [i];
        assign (arr [j], i);
        assign (item, j);
      }
    void heapify (size_t index,
                  size_t maxIndex)
      // Requires: Heap property holds for all index1 < maxIndex except parent(index1) == index
      { if (maxIndex > arr. size ())  throwError ("heapify: maxIndex");
        if (index >= maxIndex)        throwError ("heapify: index");
        for (;;)
        { size_t extr = index;
          const size_t l = left (index);
          if (   l < maxIndex
              && comp (arr [extr], arr [l]) < 0
             )
            extr = l;
          const size_t r = right (index);
          if (   r < maxIndex
              && comp (arr [extr], arr [r]) < 0
             )
            extr = r;
          if (extr == index)
            break;
          swap (index, extr);
          index = extr;
        }
      }
  public:

    // Test
    static void testStr ()    
      { const StringVector vec {"Moscow", "San Diego", "Los Angeles", "Paris"};
        Heap<const string> heap (strComp);
        for (const string& s : vec)
          heap << & s;
        while (! heap. empty ())  
        { cout << * heap. getMaximum () << endl;
          heap. deleteMaximum ();
        }
      }
  private:
    static int strComp (const void* s1,
                        const void* s2)
      { const string& s1_ = * static_cast <const string*> (s1);
        const string& s2_ = * static_cast <const string*> (s2);
        if (s1_ < s2_)  return  1;
        if (s1_ > s2_)  return -1;
        return 0;
      }
  };




///////////////////////////////////////////////////////////////////////////////

// Verbosity

bool verbose (int inc = 0);

int getVerbosity ();

inline bool getQuiet ()
  { return ! verbose (1); }


class Verbose
{
	const int verbose_old;
public:

	
	explicit Verbose (int verbose_arg);
	Verbose ();
	  // Increase verbosity
 ~Verbose ();

 
  static bool enabled ()
    { return isMainThread () /*Threads::empty ()*/; }
};



struct Unverbose 
{
	Unverbose ();
 ~Unverbose ();
};
  


struct Chronometer : Nocopy
// CPU (not astronomical) time
// Requires: no thread is used
{
  static constexpr clock_t noclock {(clock_t) -1};
  static bool enabled;
  const string name;
  clock_t time {0};
protected:
  clock_t startTime {noclock};
public:


  explicit Chronometer (const string &name_arg)
    : name (name_arg)
    {}


  bool on () const
    { return enabled && threads_max == 1; }
  bool started () const
    { return startTime != noclock; }
  void start ();
    // Requires: !started()
  void stop ();
    // Requires: started()
  void cancel ()
    { startTime = noclock; }
  void print (ostream &os) const;
};



struct Chronometer_OnePass : Nocopy
// Astronomical time
{
private:
	const string name;
  ostream &os;
  const bool addNewLine;
  const bool active;
	const time_t start;
public:
  static constexpr Color::Type color {Color::magenta};  // PAR
  
	
  explicit Chronometer_OnePass (const string &name_arg,
                                ostream &os_arg = cout,
                                bool addNewLine_arg = true,
                                bool active_arg = true);
 ~Chronometer_OnePass ();
    // Print to os
};


	

/////////////////////////////////////// cerr //////////////////////////////////////////

struct Progress : Nocopy
{
private:
	static size_t beingUsed;
public:

	const size_t n_max;
	  // 0 <=> unknown
	const bool active;
	const size_t displayPeriod;
	size_t n {0};
	string step;
	  // To report in case of throw
	

	explicit Progress (size_t n_max_arg = 0,
	                   size_t displayPeriod_arg = 1);
 ~Progress ();
    

  bool operator() (const string& step_arg = noString);
private:
	void report () const;
	  // Output: cerr
public:
  void reset ()
    { n = 0; 
      step. clear ();
    }
	static void disable ()
	  { beingUsed++; }
	static bool isUsed ()
	  { return beingUsed; }
	static bool enabled ()
	  { return ! beingUsed && ! Threads::isQuiet () && verbose (1); }
};




struct Stderr : Singleton<Stderr>
{
  bool quiet {false};

  
  template <typename T>
    Stderr& operator<< (const T& t) 
      { if (quiet)
          return *this;
        cerr << t;
        return *this;
      }
  void section (const string &title) const
    { if (quiet)
        return;
      Common_sp::section (title, false);
    }
};



struct Warning : Singleton<Warning>
{
private:
  Stderr& stderr;
  const OColor oc;
public:  

  
  explicit Warning (Stderr &stderr_arg)
    : stderr (stderr_arg)
    , oc (cerr, Color::yellow /*PAR*/, true, ! stderr. quiet)  
    { stderr << "WARNING: "; }
 ~Warning ()
    { stderr << "\n"; }
};



struct Chronometer_OnePass_cerr : Chronometer_OnePass
{
  explicit Chronometer_OnePass_cerr (const string &name_arg)
    : Chronometer_OnePass (name_arg, cerr, false, qc_on && ! getQuiet ())
    {}
};

	


////////////////////////////////////// Input ///////////////////////////////////////////

struct TextPos
// A line ends with an end-of-line (EOL)
{
  // 0-based
	int lineNum {-1};
	  // = number of lines read
private:
	static constexpr uint eol_pos {(uint) -1};
public:
  uint charNum {eol_pos};
    // In line lineNum
	
	
  bool operator== (const TextPos &other) const
    { return    lineNum == other. lineNum
    	       && charNum == other. charNum;
    }
  string str () const;

	
	void inc (bool eol_arg);
	void dec ();
	  // Requires: after inc()
	
	bool eol () const
	  { return charNum == eol_pos; }
	bool last () const
	  { return charNum + 1 == eol_pos; }


  struct Error : runtime_error
    { Error (const TextPos tp,
             const string &what,
        	   bool expected = true)
        : runtime_error ((tp. str () + what + ifS (expected, " is expected")). c_str ())
        {}
    };
};



struct Input : Root, Nocopy
{
protected:
  IFStream ifs;
  istream* is {nullptr};
    // ifs.is_open() => is = &ifs
public:
  bool eof {false};
    // End of file
protected:
	Progress prog;
public:


protected:	
  Input (const string &fName,
         uint displayPeriod)
    : ifs (fName)
    , is (& ifs)
    , prog (0, displayPeriod)  
    {}
  Input (istream &is_arg,
	       uint displayPeriod);
public:


  void reset ();
};
	


struct LineInput : Input
{
  uint lineNum {0};
	  // Number of lines read
	string line;
	  // Current line
  string commentStart;

	
	explicit LineInput (const string &fName,
          	          uint displayPeriod = 0)
    : Input (fName, displayPeriod)
    {}
  explicit LineInput (istream &is_arg,
	                    uint displayPeriod = 0)
    : Input (is_arg, displayPeriod)
    {}


	bool nextLine ();
  	// Output: line
	bool expectPrefix (const string &prefix,
	                   bool eofAllowed)
		{ if (nextLine () && trimPrefix (line, prefix))
		  	return true;  
			if (eof && eofAllowed)
				return false;
		  throw runtime_error ("No " + strQuote (prefix));
		  return false;  // dummy
		}
	string lineStr (bool add1 = true) const
	  { return "line " + to_string (lineNum + add1); }
};
	


struct PairFile : Root
{
private:
	LineInput f;
	Istringstream iss;
public:
  const bool sameAllowed;
	const bool orderNames;
  	// true => name1 < name2
	string name1;
	string name2;
	

	PairFile (const string &fName,
	          bool sameAllowed_arg,
	          bool orderNames_arg,
	          uint displayPeriod = 1000)
	  : f (fName, displayPeriod) 
	  , sameAllowed (sameAllowed_arg)
	  , orderNames (orderNames_arg)
	  {}

	  
	bool next ();
	  // Output: name1, name2
	uint getLineNum () const
	  { return f. lineNum; }
	  // Requires: after next()
};



struct CharInput : Input
{
	TextPos tp;
private:
  bool ungot {false};
public:


	explicit CharInput (const string &fName,
              	      uint displayPeriod = 0)
    : Input (fName, displayPeriod)
    {}
  explicit CharInput (istream &is_arg,
	                    uint displayPeriod = 0)
    : Input (is_arg, displayPeriod)
    {}


	char get ();
	  // Update: tp
	void unget ();
	  // Requires: To be followed by get()
  string getLine ();
    // Postcondition: tp.eol()
	  
  [[noreturn]] void error (const string &what,
	                         bool expected = true) const
		{ throw TextPos::Error (tp, what, expected); }
};
	


struct Token : Root
{
	enum Type { eName
	          , eDelimiter  
	          , eText
	          , eInteger   // 10-based or 16-based: 0xNNNN...
	          , eDouble
	          , eDateTime  // Example: 2018-08-13T16:12:54.487
	          };
  // Valid if !empty()
	TextPos tp;
    // = CharInput::tp
	Type type {eDelimiter};
	bool dashInName {false};
	string name;
	  // eName => !empty(), printable()
	  // eText => embracing quote's are removed
	  // eDilimiter => name.size() = 1
	char quote {'\0'};
	  // '\0': no quote
	long long n {0};
	double d {0.0};
  // Valid if eDouble
	streamsize decimals {0};
	bool scientific {false};	 

	  
	Token () = default;
	Token (const string& name_arg,
	       bool dashInName_arg)
	  : type (eName)
	  , dashInName (dashInName_arg)
	  , name (name_arg)
	  {}
	Token (const string& name_arg,
	       char quote_arg)
	  : type (eText)
	  , name (name_arg)
	  , quote (quote_arg)
	  {}
	explicit Token (long long n_arg)
	  : type (eInteger)
	  , n (n_arg)
	  {}
	explicit Token (double d_arg)
	  : type (eDouble)
	  , d (d_arg)
	  {}
	explicit Token (char delimiter_arg)
	  : type (eDelimiter)
	  , name (1, delimiter_arg)
	  {}
	Token (CharInput &in,
	       bool dashInName_arg,
	       bool consecutiveQuotesInText,
	       bool negativeInteger)
	  { readInput (in, dashInName_arg, consecutiveQuotesInText, negativeInteger); }
	Token (CharInput &in,
	       Type expected,
	       bool dashInName_arg,
	       bool consecutiveQuotesInText,
	       bool negativeInteger)
    { readInput (in, dashInName_arg, consecutiveQuotesInText, negativeInteger);
    	if (empty ())
 			  in. error ("No token", false); 
    	if (type != expected)
 			  in. error (type2str (type)); 
    }
private:  
	void readInput (CharInput &in,
	                bool dashInName_arg,
	                bool consecutiveQuotesInText,
	                bool negativeInteger);  
	  // Input: consecutiveQuotesInText means that '' = '
    // Update: in: in.charNum = last character of *this
    // --> TokenInput ??
public:
	void qc () const override;
	void saveText (ostream &os) const override;
	bool empty () const override
	  { return type == eDelimiter && name. empty (); }
	void clear () override
	  { *this = Token (); }


  [[noreturn]] void error (const string &what,
                           bool expected = true) const
		{ throw TextPos::Error (tp, what, expected); }  
	static string type2str (Type type) 
	  { switch (type)
	  	{ case eName:      return "name";
	  		case eText:      return "text";
	  		case eInteger:   return "integer";
	  		case eDouble:    return "double";
	  		case eDelimiter: return "delimiter";
	  		case eDateTime:  return "datetime";
	  	}
	  	throwf ("Token::type2str()");
	  	return noString;
	  }
	static Type str2type (const string &s)
	  { if (s == "name")      return eName;
	    if (s == "text")      return eText;
	    if (s == "integer")   return eInteger;
	    if (s == "double")    return eDouble;
	    if (s == "delimiter") return eDelimiter;
	    if (s == "datetime")  return eDateTime;
	    throwf ("Token::str2type()");
	    return eName;  // dummy
	  }      
	bool isName (const string &s) const
	  { return ! empty () && type == eName && name == s; }
	bool isNameText (const string &s) const
	  { return ! empty () && (type == eName || type == eText) && name == s; }
	bool isInteger (long long n_arg) const
	  { return ! empty () && type == eInteger && n == n_arg; }
	bool isDouble (double d_arg) const
	  { return ! empty () && type == eDouble && d == d_arg; }
	bool isDelimiter (char c) const
	  { return ! empty () && type == eDelimiter && name [0] == c; }
	  
	bool operator== (const Token &other) const
	  { if (type != other. type)
	      return false;
	    switch (type)
	    { case eName:
	      case eText:
	      case eDateTime:
	      case eDelimiter: return name == other. name;
	      case eInteger:   return n    == other. n;
	      case eDouble:    return d    == other. d;
	    }
	    return false;
	  }
	bool operator< (const Token &other) const;
	
	void toNumberDate ();
	  // Try to convert eName or eText to: eInteger, eDouble, eDateTime
};



struct TokenInput : Root
{
private:
  CharInput ci;
public:
  const char commentStart;
  const bool dashInName;
  const bool consecutiveQuotesInText;
    // Two quotes encode one quote
  const bool negativeInteger;
private:
  Token last;
	TextPos tp;
public:


  explicit TokenInput (const string &fName,
                       char commentStart_arg = '\0',
                       bool dashInName_arg = false,
                       bool consecutiveQuotesInText_arg = false,
                       bool negativeInteger_arg = false,
                       uint displayPeriod = 0)
    : ci (fName, displayPeriod)
    , commentStart (commentStart_arg)
    , dashInName (dashInName_arg)
    , consecutiveQuotesInText (consecutiveQuotesInText_arg)
    , negativeInteger (negativeInteger_arg)
    {}
  explicit TokenInput (istream &is_arg,
                       char commentStart_arg = '\0',
                       bool dashInName_arg = false,
                       bool consecutiveQuotesInText_arg = false,
                       bool negativeInteger_arg = false,
                       uint displayPeriod = 0)
    : ci (is_arg, displayPeriod)
    , commentStart (commentStart_arg)
    , dashInName (dashInName_arg)
    , consecutiveQuotesInText (consecutiveQuotesInText_arg)
    , negativeInteger (negativeInteger_arg)
    {}


  [[noreturn]] void error (const string &what,
	                         bool expected = true) const
		{ throw TextPos::Error (tp, what, expected); }  

  Token get ();
    // Return: empty() <=> EOF
  Token getXmlText ();
    // Last character read is '<', must be followed by '/'
  Token getXmlComment ();
    // -- ... -->
  Token getXmlProcessingInstruction ();
    // ... &>
  Token getXmlMarkupDeclaration ();
    // ... >
  char getNextChar (bool unget);
    // Return: '\0' <=> EOF
    // Invokes: if (unget) ci.unget()
	void get (const string &expected)
    { const Token t (get ());
    	if (! t. isNameText (expected))
   			t. error (Token::type2str (Token::eName) + " " + strQuote (expected)); 
    }
	void get (int expected)
    { const Token t (get ());
    	if (! t. isInteger (expected))
  			t. error (Token::type2str (Token::eInteger) + " " + to_string (expected)); 
    }
	void get (double expected)
    { const Token t (get ());
    	if (! t. isDouble (expected))
   			t. error (Token::type2str (Token::eDouble) + " " + toString (expected)); 
    }
	void get (char expected)
    { const Token t (get ());
      if (! t. isDelimiter (expected))
   			error (Token::type2str (Token::eDelimiter) + " " + strQuote (toString (expected), '\'')); 
    }
  void setLast (Token &&t)
    { last = std::move (t); 
      tp = last. tp;
    }
  bool getNext (char expected)
    { Token token (get ());
      if (! token. isDelimiter (expected))
      { setLast (std::move (token));
      	return false;
      }
      return true;
    }
};




struct BraceInput : TokenInput
{
  static constexpr char commentC {'#'};
  static constexpr const char* commentS {"comment"};
  static constexpr char endChar {';'};


  explicit BraceInput (const string &fName)
    : TokenInput (fName, commentC)
    {}
  explicit BraceInput (istream &is_arg)
    : TokenInput (is_arg, commentC)
    {}


  static string endS ()
    { return string (1, endChar); }
  void skipComment ();
};




///////////////////////////////////// Json //////////////////////////////////////////
		
extern unique_ptr<JsonMap> jRoot;



struct Json : Root, Nocopy  // Heaponly
{
protected:
  Json (JsonContainer* parent,
        const string& name);
  Json () = default;
public:  
  void saveText (ostream& os) const override = 0;
  

  virtual const JsonNull* asJsonNull () const
    { return nullptr; }  
  virtual const JsonString* asJsonString () const
    { return nullptr; }  
  virtual const JsonInt* asJsonInt () const
    { return nullptr; }  
  virtual const JsonDouble* asJsonDouble () const
    { return nullptr; }  
  virtual const JsonBoolean* asJsonBoolean () const
    { return nullptr; }  
  virtual const JsonArray* asJsonArray () const
    { return nullptr; }  
  virtual const JsonMap* asJsonMap () const
    { return nullptr; }  

protected:
  static string toStr (const string& s);
  static void parse (CharInput &in,
                     const Token& firstToken,
                     JsonContainer* parent,
                     const string& name);
public:    
  string getString () const;
    // Requires: JsonString
  long long getInt () const;
    // Requires: JsonInt
  double getDouble () const;
    // Requires: JsonDouble
  bool getBoolean () const;
    // Requires: JsonBoolean
  const Json* at (const string& name_arg) const;
    // Requires: JsonMap
  const Json* at (size_t index) const;
    // Requires: JsonArray
  size_t getSize () const;
    // Requires: JsonArray
};



struct JsonNull final : Json
{
  explicit JsonNull (JsonContainer* parent,
                     const string& name = noString)
    : Json (parent, name)
    {}    
  void saveText (ostream& os) const final
    { os << "null"; }


  const JsonNull* asJsonNull () const final
    { return this; }  
};



struct JsonString final : Json
{
  const string s;


  JsonString (const string& s_arg,
              JsonContainer* parent,
              const string& name = noString)
    : Json (parent, name)
    , s (s_arg)
    {}
  void saveText (ostream& os) const final
    { os << toStr (s); }


  const JsonString* asJsonString () const final
    { return this; }  
};



struct JsonInt final : Json
{
  const long long n;
  
  JsonInt (long long n_arg,
           JsonContainer* parent,
           const string& name = noString)
    : Json (parent, name)
    , n (n_arg)
    {}
  void saveText (ostream& os) const final
    { os << n; }


  const JsonInt* asJsonInt () const final
    { return this; }  
};



struct JsonDouble final : Json
{
  const double n;
  const streamsize decimals;
  bool scientific {false};


  JsonDouble (double n_arg,
              streamsize decimals_arg,
              JsonContainer* parent,
              const string& name = noString)
    : Json (parent, name)
    , n (n_arg)
    , decimals (decimals_arg == numeric_limits<streamsize>::max() ? double2decimals (n_arg) : decimals_arg)
    {}
    // decimals_arg = -1: default
  void saveText (ostream& os) const final
    { const ONumber on (os, (streamsize) decimals, scientific);
    	if (isfinite (n))
        os << n; 
      else
        os << "null";  // NaN
    }      


  const JsonDouble* asJsonDouble () const final
    { return this; }  
};



struct JsonBoolean final : Json
{
  const bool b;
  

  JsonBoolean (bool b_arg,
               JsonContainer* parent,
               const string& name = noString)
    : Json (parent, name)
    , b (b_arg)
    {}
  void saveText (ostream& os) const final
    { os << (b ? "true" : "false"); }


  const JsonBoolean* asJsonBoolean () const final
    { return this; }  
};



struct JsonContainer : Json
{
protected:
  JsonContainer (JsonContainer* parent,
                 const string& name)
    : Json (parent, name)
    {}
  JsonContainer () = default;
public:  
};



struct JsonArray final : JsonContainer
{
  friend struct Json;
private:
  VectorOwn<Json> data;
public:
  

  explicit JsonArray (JsonContainer* parent,
                      const string& name = noString)
    : JsonContainer (parent, name)
    {}
private:
  JsonArray (CharInput& in,
             JsonContainer* parent,
             const string& name);
public:
  void saveText (ostream& os) const final;


  const JsonArray* asJsonArray () const final
    { return this; }
    
  size_t size () const
    { return data. size (); }
};



struct JsonMap final : JsonContainer
{
  friend struct Json;
private:
  typedef  map <string, const Json*>  Map;
  Map data;
public:
  
  
  explicit JsonMap (JsonContainer* parent,
                    const string& name = noString)
    : JsonContainer (parent, name)
    {}
  JsonMap ();
    // Output: jRoot = this
  explicit JsonMap (const string &fName);
private:
  JsonMap (CharInput& in,
           JsonContainer* parent,
           const string& name)
    : JsonContainer (parent, name)
    { parse (in); }
  void parse (CharInput& in);
public:
 ~JsonMap ();
  void saveText (ostream& os) const final;


  const JsonMap* asJsonMap () const final
    { return this; }
    
  StringVector getKeys () const
    { StringVector keys;  keys. reserve (data. size ());
      for (const auto& it : data)
        keys << it. first;
      return keys;
    }
};



template <typename T>
  Json* vec2json (const Vector<T> &vec,
                  JsonContainer* parent,
                  const string& name = noString) 
    { auto j = new JsonArray (parent, name);
      for (const T& t : vec)
        t. toJson (j);
      return j;
    }



///////////////////////////////////////////////////////////////////////////

struct ItemGenerator
{
  Progress prog;
  
  
protected:
  ItemGenerator (size_t progress_n_max,
	               size_t progress_displayPeriod)
	  : prog (progress_n_max, progress_displayPeriod)
	  {}
public:
  virtual ~ItemGenerator ()
    {}
  
  
  virtual bool next (string &item) = 0;
    // Return: false <=> end of items
    // Output: item; may be empty()
    // Invokes: prog()
};



struct FileItemGenerator final : ItemGenerator, Nocopy
{
private:
  string fName;
  IFStream f;
public:
  const bool tsv;
  
  
  FileItemGenerator (size_t progress_displayPeriod,
                     const string& fName_arg,
                     bool tsv_arg);
  
  
  bool next (string &item) final;
};

  

#ifndef _MSC_VER
struct RawDirItemGenerator final : ItemGenerator, Nocopy
{
private:
  string dirName;
  struct Imp;  
  Imp* imp {nullptr};
  unique_ptr<RawDirItemGenerator> dig;
public:
  const bool large;
  
  
  RawDirItemGenerator (size_t progress_displayPeriod,
                       const string& dirName_arg,
                       bool large_arg);
 ~RawDirItemGenerator ();

  
  bool next (string &item) final;
private:
  bool next_ (string &item,
              bool report);
public:
  StringVector toVector ();
};



struct DirItemGenerator final : ItemGenerator
{
private:
  StringVector vec;
  size_t index {0};
public:
  
  
  DirItemGenerator (size_t progress_displayPeriod,
                    const string& dirName,
                    bool large)
    : ItemGenerator (0, progress_displayPeriod)
    , vec (RawDirItemGenerator (0, dirName, large). toVector ())
    {}

  
  bool next (string &item) final
    { if (index == vec. size ())
        return false;
      item = vec [index];
      prog (item);
      index++;      
      return true;
    }
};
#endif

  

struct NumberItemGenerator final : ItemGenerator
{
private:
  size_t index {0};
public:
  
  
  NumberItemGenerator (size_t progress_displayPeriod,
                       size_t n)
    : ItemGenerator (n, progress_displayPeriod)
    {}
  
  
  bool next (string &item) final
    { if (index == prog. n_max)
        return false;
      index++;
      item = to_string (index);
      prog ();
      return true;
    }
};

  


////////////////////////////////// Application //////////////////////////////////////////

struct SoftwareVersion : Root
{
  uint major {0};  // there is ::major()
  uint minor {0};  // there is ::minor()
  uint patch {0};
  

  SoftwareVersion () = default;
  SoftwareVersion (uint major_arg,
                   uint minor_arg,
                   uint patch_arg)
    { major = major_arg;
      minor = minor_arg;
      patch = patch_arg;
    } 
  explicit SoftwareVersion (const string &fName);
  explicit SoftwareVersion (istream &is,
                            bool minorOnly = false);
private:
  void init (string &&s,
             bool minorOnly);
public:
  void saveText (ostream &os) const final
    { os << major << '.' << minor << '.' << patch; }   
    
    
  bool operator< (const SoftwareVersion &other) const;
  bool operator== (const SoftwareVersion &other) const
    { return    major == other. major
             && minor == other. minor
             && patch == other. patch;
    }
  bool operator<= (const SoftwareVersion &other) const
    { return operator< (other) || operator== (other); }
    
  string getMinor () const
    { return to_string (major) + "." + to_string (minor); }
};
  
  

struct DataVersion : Root
{
  uint year  {0};
  uint month {0};
  uint day   {0};
  uint num   {0};
  

  DataVersion () = default;
  DataVersion (uint year_arg,
               uint month_arg,
               uint day_arg,
               uint num_arg)
    : year  (year_arg)
    , month (month_arg)
    , day   (day_arg)
    , num   (num_arg)
    {} 
  explicit DataVersion (const string &fName);
  explicit DataVersion (istream &is);
private:
  void init (string &&s);
public:
  void saveText (ostream &os) const final
    { os << year 
         << '-' << std::setfill ('0') << std::setw (2) << month 
         << '-' << std::setfill ('0') << std::setw (2) << day 
         << '.' << num; 
    }   
    
    
  bool operator< (const DataVersion &other) const;
  bool operator== (const DataVersion &other) const
    { return    year  == other. year
             && month == other. month
             && day   == other. day
             && num   == other. num;
    }
  bool operator<= (const DataVersion &other) const
    { return operator< (other) || operator== (other); }
};
  
  

struct Application : Singleton<Application>, Root
// Usage: int main (argc, argv) { ThisApplication /*:Application*/ app; return app. run (argc, argv); }
{  
  const string description;
  string version {"0.0.0"};
  // Without "https://"
  string documentationUrl;  
  string updatesUrl;  
  //
  string updatesDoc;
  const bool needsArg;
  const bool gnu;
  const bool threadsUsed;
  string execDir;
    // Ends with '/'
    // Physically real directory of the software
  static constexpr const char* helpS {"help"};
  static constexpr const char* versionS {"version"};
  
protected:
  struct Positional;  // forward
  struct Key;         // forward
  struct Arg : Named
  {
    const string description;
      // !empty()
    string value;
  protected:
    Arg (const string &name_arg,
         const string &description_arg);
  public:
  	void qc () const override;
    virtual const Positional* asPositional () const
      { return nullptr; }
    virtual const Key* asKey () const
      { return nullptr; }
  };
  struct Positional final : Arg
  {
    Positional (const string &name_arg, 
                const string &description_arg)
      : Arg (name_arg, description_arg)
      {}
    void saveText (ostream &os) const override
      { os << '<' << name << '>'; }
    const Positional* asPositional () const final
      { return this; }
  };  
  struct Key final : Arg
  {
  	const Application& app;
    const bool flag;
    string requiredGroup;
  	const char acronym;
  	  // '\0' <=> no acronym
    const string var;
      // For help
    string defaultValue;
    Key (const Application &app_arg,
         const string &name_arg,
         const string &description_arg,
         const string &defaultValue_arg,
         char acronym_arg,
         const string &var_arg)
      : Arg (name_arg, description_arg)
      , app (app_arg)
      , flag (false)
      , acronym (acronym_arg)
      , var (var_arg)
      , defaultValue (defaultValue_arg)
      {}
      // "$BASE" in defaultValue means `dirname $0`
    Key (const Application &app_arg,
         const string &name_arg,
         const string &description_arg,
         char acronym_arg)
      : Arg (name_arg, description_arg)
      , app (app_arg)
      , flag (true)
      , acronym (acronym_arg)
      {}
  	void qc () const override;
    void saveText (ostream &os) const override;
    string getShortHelp () const;
    const Key* asKey () const final
      { return this; }
  };
  List<Positional> positionalArgs;
  List<Key> keyArgs;
  map<string/*Arg::name*/,const Arg*> name2arg;
  map<char/*Arg::name[0]*/,const Key*> char2arg;
    // Valid if gnu
public:
  

protected:
  explicit Application (const string &description_arg,
                        bool needsArg_arg = true,
                        bool gnu_arg = false,
                        bool threadsUsed_arg = false)   
    : description (description_arg)
    , needsArg (needsArg_arg)
    , gnu (gnu_arg)
    , threadsUsed (threadsUsed_arg)
    {}
    // To invoke: addKey(), addFlag(), addPositional(), setRequiredGroup()
  // <Command-line parameters> ::= <arg>*
  //   <arg> ::= <positional> | <key> | <flag>
  //   <positional> ::= <string>
  //   <key> ::= -<name> <value> | -<name>=<value>
  //   <flag> ::= -<name>
  // acronym = '\0' <=> no acronym
  void addKey (const string &name, 
               const string &argDescription,
               const string &defaultValue = noString,
               char acronym = '\0',
               const string &var = noString);
    // [-<name> <defaultValue>]
    // gnu: [--<name> <var>]
  void addFlag (const string &name,
                const string &argDescription,
                char acronym = '\0');
    // [-<name>]
  void addPositional (const string &name,
                      const string &argDescription);
    // Requires: gnu => name is uppercased
  void setRequiredGroup (const string &keyName,
                         const string &requiredGroup);
private:
	void addDefaultArgs ();
	void qc () const final;
	Key* getKey (const string &keyName) const;
	  // Return: !nullptr
	void setPositional (List<Positional>::iterator &posIt,
	                    const string &value);
public:
 ~Application ()
   { if (logPtr)
  	 { delete logPtr;
  	   logPtr = nullptr;
     }
   }


protected:
  string getArg (const string &name) const;
    // Input: keyArgs, where Key::flag = false, and positionalArgs
  uint arg2uint (const string &name) const
    { uint n = 0;
    	if (! str2<uint> (getArg (name), n))
    	  throw runtime_error ("Cannot convert " + strQuote (name) + " to non-negative number"); 
    	return n;
    }
  double arg2double (const string &name) const
    { double d = 0.0;
    	if (! str2<double> (getArg (name), d))
    	  throw runtime_error ("Cannot convert " + strQuote (name) + " to real number"); 
    	return d;
    }
  bool getFlag (const string &name) const;
    // Input: keyArgs, where Key::flag = true
  string key2shortHelp (const string &name) const;
  string getProgramDirName () const
    { return getDirName (programArgs. front ()); }
protected:
  virtual void initEnvironment ();
  virtual void initVar ()
    {}
  string getInstruction (bool screen) const;
  string getDocumentation (bool screen) const;
  string getUpdates (bool screen) const;
  virtual string getHelp (bool screen) const;
    // Requires: must be printed to cout
  string makeKey (const string &param,
                  const string &value) const
    { if (value. empty ())
        return noString;
      return "  -" + ifS (gnu, "-") + param + " " + shellQuote (value); 
    }
public:
  int run (int argc, 
           const char* argv []);
    // Invokes: body()
    // if runtime_error then exit(1) else errorExit()
private:
  virtual void body () const = 0;
    // Invokes: initEnvironment()
};



#ifndef _MSC_VER
struct ShellApplication : Application
// Requires: SHELL=bash
{
protected:
  // Environment
  const bool useTmp;
  string tmp;
    // Temporary directory: ($TMPDIR or "/tmp") + "/" + programName + "XXXXXX"
    // If log is used then tmp is printed in the log file and the temporary files are not deleted 
    // !empty() => useTmp
    // See https://systemd.io/TEMPORARY_DIRECTORIES/
  mutable KeyValue prog2dir;
  mutable Stderr stderr;
  time_t startTime {0};
public:
  
  
  ShellApplication (const string &description_arg,
                    bool needsArg_arg,
                    bool threadsUsed_arg,
                    bool gnu_arg,
                    bool useTmp_arg)
    : Application (description_arg, needsArg_arg, gnu_arg, threadsUsed_arg)
    , useTmp (useTmp_arg)
    {}
 ~ShellApplication ();
   // Invokes: removeDirectory(tmp) if !logPtr


protected:
  void initEnvironment () override;
  void initVar () override;
    // stderr << "Running: " << getCommandLine ()
    // threads_max correction
  string getHelp (bool screen) const override;
private:
  void body () const final;
  virtual void shellBody () const = 0;
protected:
  static bool emptyArg (const string &s)
    {	return s. empty () || s == "\'\'"; }
  void findProg (const string &progName) const;
    // Output: prog2dir
  string fullProg (const string &progName) const;
    // Return: shellQuote (directory + progName) + ' '
    // Requires: After findProg(progName)
  string exec2str (const string &cmd,
                   const string &tmpName,
                   const string &logFName = noString) const;
    // Return: $( cmd > <tmp>/tmpName && cat <tmp>/tmpName )
    // Requires: cmd produces one line
  string uncompress (const string &quotedFName,
                     const string &suffix) const;
    // Return: quotedFName or tmp + "/" + suffix
    // Invokes: exec("gunzip")
  string getBlastThreadsParam (const string &blast,
                               size_t threads_max_max) const;
    // Input: blast: "blastp", "blastx", etc.
    // Invokes: fullProg(blast)
};
#endif



}



#endif

