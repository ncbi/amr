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
  #pragma warning (disable : 4290)  // C++ exception specification ignored except to indicate a function is not __declspec(nothrow)
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

#include <ctime>
#include <cstring>
#include <string>
#include <stdexcept>
#include <limits>
#include <array>
#include <list>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <algorithm>

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


typedef  unsigned char  uchar; 
typedef  unsigned int   uint; 
typedef  unsigned long  ulong; 



bool initCommon ();
  // Module initialization
  // Invoked automaticallly

extern vector<string> programArgs;
extern string programName;

string getCommandLine ();

extern ostream* logPtr;

extern bool qc_on;
extern ulong seed_global;
  // >= 1

// thread
extern size_t threads_max;
  // >= 1
extern thread::id main_thread_id;
bool isMainThread ();


void errorExit (const char* msg,
                bool segmFault = false);
  // Input: programArgs, programName, logptr
	// Update: *logPtr
	// Invokes: if segmFault then abort() else exit(1)

inline void errorExitStr (const string &msg)
  { errorExit (msg. c_str ()); }



struct Nocopy
{
protected:
	Nocopy () = default;
  Nocopy (const Nocopy &) = delete;
  Nocopy (Nocopy &&) = delete;
  Nocopy& operator= (const Nocopy &) = delete;
};



class ONumber
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



struct Chronometer : Nocopy
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
  void start ()
    { if (! on ())
        return;
      if (startTime != noclock)
        throw logic_error ("Chronometer \""  + name + "\" is not stopped");
      startTime = clock (); 
    }  
  void stop () 
    { if (! on ())
        return;
      if (startTime == noclock)
        throw logic_error ("Chronometer \"" + name + "\" is not started");
      time += clock () - startTime; 
      startTime = noclock;
    }

  void print (ostream &os) const
    { if (! on ())
        return;       
      os << "CHRON: " << name << ": ";
      const ONumber onm (os, 2, false);
      os << (double) time / CLOCKS_PER_SEC << " sec." << endl;
    }
};



struct Chronometer_OnePass : Nocopy
{
	const string name;
	const time_t start;
	
  explicit Chronometer_OnePass (const string &name_arg)
    : name (name_arg)
    , start (time (nullptr))
    {}
 ~Chronometer_OnePass ()
    { if (uncaught_exception ())
        return;
      const time_t stop = time (nullptr);
      cout << "CHRON: " << name << ": ";
      const ONumber onm (cout, 0, false);
      cout << difftime (stop, start) << " sec." << endl;
      cout << endl;
    }
};
	
	

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
      throw runtime_error ("Dereferencing nullptr");
    }



// bool

inline void toggle (bool &b)
  { b = ! b; }

inline int getSign (bool b)
  { return b ? 1 : -1; }

inline bool boolPow (bool x, bool power)
  { return power ? x : ! x; }



// ebool - extended bool

enum ebool {EFALSE = false, 
            ETRUE = true, 
            UBOOL = true + 1};

inline ebool toEbool (bool b)
  { return b ? ETRUE : EFALSE; }


inline bool operator<= (ebool a, ebool b)
  { static constexpr char rank [3/*ebool*/] = {0, 2, 1};
  	return rank [a] <= rank [b];
  }

inline void toggle (ebool &b)
  { if (b == ETRUE)
      b = EFALSE;
    else if (b == EFALSE)
      b = ETRUE;
  }



inline void advance (size_t &index, 
                     size_t size)
  // Update: index: < size
  { index++; if (index == size) index = 0; }
  
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
  
template <typename T/*:integer*/> 
  inline bool even (T x)
    { static_assert (numeric_limits<T>::is_integer, "Must be integer");
      return x % 2 == 0; 
    }

inline bool divisible (uint n,
                       uint divisor)
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

uint powInt (uint a,
             uint b);
  // Return: a^b
  // Time: O(log(b))



const size_t NO_INDEX = SIZE_MAX;



// char

inline bool isAlpha (char c)
  { return strchr ("abcdefghijklmnopqrstuvwxyz", tolower (c)); }

inline bool isDigit (char c)
  { return strchr ("0123456789", c); }
  
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
 


/* Usage:
  for (Iter<T> iter (t); iter. next (); )
    ...
*/
template <typename T>
struct Iter : Nocopy
{
private:
  T& t;
  typename T::iterator itNext;
public:
  typename T::iterator it;
    

  Iter (T &t_arg)
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
    { typename T::value_type val = *it;
      itNext = t. erase (it); 
      return val;
    }
  void insert (const typename T::value_type &val)
    { itNext = t. insert (itNext, val) + 1; }
    // Opposite to erase()
};



template <typename T, typename U>
  inline ostream& operator<< (ostream &os,
                              const pair<T,U> &p) 
    { os << p. first << '\t' << p. second;
      return os;
    }



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
  Pair () = default;
    
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
};



// STL algorithms

template <typename To, typename From>
  inline void insertAll (To &to,
                         const From &from)
    { to. insert (to. begin (), from. begin (), from. end ()); }

template <typename To, typename From>
  inline void insertIter (To &to,
                          const From &from)
    { for (const auto x : from)
        to << x;
    }

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
    { for (const auto x : t)
        if (u. find (static_cast <typename U::value_type> (x)) != u. end ())
          return true;
      return false;
    }

template <typename Key, typename Value, typename KeyParent>
  inline bool contains (const map <Key, Value> &m,
                        const KeyParent& key)
    { return m. find (key) != m. end (); }

template <typename Key, typename Value, typename KeyParent>
  inline bool contains (const unordered_map <Key, Value> &m,
                        const KeyParent& key)
    { return m. find (key) != m. end (); }

template <typename Key, typename Value, typename KeyParent>
  bool find (const map <Key, Value> &m,
             const KeyParent& key,
             Value &value)
    // Return: success
    // Output: value, if Return
    { const auto it = m. find (key);
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
    { const auto it = m. find (key);
    	if (it == m. end ())
    		return false;
    	value = it->second; 
    	return true;
    }

template <typename Key, typename Value, typename KeyParent>
  const Value* findPtr (const map <Key, Value> &m,
                        const KeyParent& key)
    { const auto it = m. find (key);
    	if (it == m. end ())
    		return nullptr;
    	return & it->second; 
    }

template <typename Key, typename Value, typename KeyParent>
  const Value* findPtr (const unordered_map <Key, Value> &m,
                        const KeyParent& key)
    { const auto it = m. find (key);
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
  inline long count_if (T &t, UnaryPredicate pred)
    { return std::count_if (t. begin (), t. end (), pred); }




template <typename T>
  ostream& operator<< (ostream &os,
                       const list<T> &ts)
    { for (const auto& t : ts)
        os << t << endl;
      return os;
    }



struct Istringstream : istringstream
{
  Istringstream () = default;
  void reset (const string &s)
    { clear ();
      str (s);
    }
};



template <typename T>
struct List : list<T>
{
private:
	typedef  list<T>  P;
public:


	List () = default;
  template <typename U/*:<T>*/>
    explicit List (const list<U> &other)
      { *this << other; }
  template <typename U/*:<T>*/>
    explicit List (const vector<U> &other)
      { *this << other; }

	  
	T at (size_t index) const
	  { size_t i = 0;
	  	for (const T& t : *this)
	  	  if (i == index)
	  	  	return t;
	  	  else
	  	  	i++;
	  	throw runtime_error ("List index is out of range");
	  }
	size_t find (const T &t) const
	  { size_t i = 0;
	  	for (const T& item : *this)
	  	  if (item == t)
	  	  	return i;
	  	  else
	  	  	i++;
	  	return NO_INDEX;
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
  List<T>& operator<< (T t) 
    { P::push_back (t); 
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
    { const T t = P::front ();
      P::pop_front ();
      return t;
    }
  T popBack ()
    { const T t = P::back ();
      P::pop_back ();
      return t;
    }
};



// string

extern const string noString;

inline string ifS (bool cond,
                   const string &s)
  { return cond ? s : string (); }

inline string nvl (const string& s,
                   const string& nullS = "-")
  { return s. empty () ? nullS : s; }
  	
inline bool isQuoted (const string &s,
                      char quote = '\"')
  { return ! s. empty () && s [0] == quote && s [s. size () - 1] == quote; }

string strQuote (const string &s,
                 char quote = '\"');

inline string unQuote (const string &s)
  { return s. substr (1, s. size () - 2); }

inline string prepend (const string &prefix,
                    	 const string &s)
  { if (s. empty ())
  	  return string ();
  	return prefix + s;
  }

inline bool isSpace (char c)
  { return c > '\0' && c <= ' ' && isspace (c); }

bool strBlank (const string &s);

template <typename T>
  string toString (const T t)
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
         )
        throw runtime_error ("Converting " + strQuote (s));
      return i;
    }

template <typename T>
  bool str2 (const string &s,
             T &t)
    { try { t = str2<T> (s); return true; } 
        catch (...) { return false; } 
    }

inline bool isLeft (const string &s,
                    const string &left)
  { return s. substr (0, left. size ()) == left; }

bool isRight (const string &s,
              const string &right);

bool trimPrefix (string &s,
                 const string &prefix);
  // Return: success

bool trimSuffix (string &s,
                 const string &suffix);
  // Return: success

void trimSuffixNonAlphaNum (string &s);

bool trimTailAt (string &s,
                 const string &tailStart);
  // Return: trimmed
  // Update: s

bool goodName (const string &name);

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

string::const_iterator stringInSet (const string &s,
                    	           	  const string &charSet);
  // Return: != s.end() => *Return is not in charSet

size_t strCountSet (const string &s,
		                const string &charSet);

void strDeleteSet (string &s,
		               const string &charSet);

void trimLeading (string &s);

void trimTrailing (string &s);

inline void trim (string &s)
  { trimTrailing (s);
    trimLeading (s); 
  }

void trimLeading (string &s,
                  char c);

void trimTrailing (string &s,
                   char c);

inline void trim (string &s,
                  char c)
  { trimTrailing (s, c);
    trimLeading  (s, c); 
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

void replace (string &s,
              char from,
              char to);
  
void replace (string &s,
              const string &fromChars,
              char to);
  
void replaceStr (string &s,
                 const string &from,
                 const string &to);
  // Replaces "from" by "to" in s from left to right
  // The replacing "to" is skipped if it contains "from"
  // Requires: !from.empty()

string to_c (const string &s);
  // " --> \", etc.

void collapseSpace (string &s);
  
string str2streamWord (const string &s,
                       size_t wordNum);
  // Return: May be string()
  // Input: wordNum: 0-based  

string str2sql (const string &s);
  // Return: ' s ' (replacing ' by '')

string sql2escaped (const string &s);


string findSplit (string &s,
                  char c = ' ');
	// Return: prefix of s+c before c
	// Update: s

string rfindSplit (string &s,
                   char c = ' ');
	// Return: suffix of c+s after c
	// Update: s

void reverse (string &s);

List<string> str2list (const string &s,
                       char c = ' ');
  // Invokes: findSplit(s,c)

string list2str (const List<string> &strList,
                 const string &sep = " ");


const char fileSlash = 
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
  		return string ();
  	return path. substr (0, pos + 1);
  }

inline bool isDirName (const string &path)
  { return isRight (path, "/"); }

bool fileExists (const string &fName);

#ifndef _MSC_VER
  bool directoryExists (const string &dirName);
#endif

streampos getFileSize (const string &fName);

size_t strMonth2num (const string& month);
// Input: month: "Jan", "Feb", ... (3 characters)


// istream

bool getChar (istream &is,
              char &c);
  // Output: s if (bool)Return

void skipLine (istream &is);

void readLine (istream &is,
               string &s);
  // Output: s

string getToken (istream &is,
                 const string &skip,
                 const string &delimeters);
  // Return: empty() <=> eof

inline void pressAnyKey ()
  { cout << "Press any key..."; char c; cin >> c; }



inline streamsize double2decimals (double r)
  { return r ? (streamsize) max<long> (0, (long) (ceil (- log10 (fabs (r)) + 1))) : 0; }



struct Rand
// Numerical Recipes in C, p. 279
{
	static const long /*actually ulong*/ max_;
private:
	long seed;
	  // 0 < seed < max_
public:
	
	explicit Rand (ulong seed_arg = 1)
  	{ setSeed (seed_arg); }
	void setSeed (ulong seed_arg)
    {	seed = (long) (seed_arg % (ulong) max_);
	    qc ();
    }
	  
	ulong get (ulong max);
    // Return: in 0 .. max - 1
    // Input: 0 < max <= max_
  double getProb ()
    { return (double) get ((ulong) max_) / ((double) max_ - 1); }
    // Return: [0,1]
private:
	void qc () const;
	void run ();	
};



template <typename T>
struct Singleton : Nocopy
{
private:
	static bool beingRun;
protected:
	Singleton ()
	  { if (beingRun)
	  	  throw runtime_error ("Singleton");
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
      { *ptr = t; }
  };



//

bool verbose (int inc = 0);



// Threads

struct Threads : Singleton<Threads>
// Usage: { Threads th; th << ...; main_thread_process(); }
{
private:
	static size_t threadsToStart;
	vector<thread> threads;
	const bool quiet;
public:

  explicit Threads (size_t threadsToStart_arg,
                    bool quiet_arg = false);
 ~Threads ();
  	
	static bool empty () 
	  { return ! threadsToStart; }
	Threads& operator<< (thread &&t)
	  { if (threads. size () >= threadsToStart)
	  	  throw logic_error ("Too many threads created");
	  	try { threads. push_back (move (t)); }
	  	  catch (const exception &e) 
	  	    { throw runtime_error (string ("Cannot start thread\n") + e. what ()); }
	  	return *this;
	  }
	size_t getAvailable () const
	  { return threadsToStart - threads. size (); }
};



template <typename Func, typename Res, typename... Args>
  void arrayThreads (const Func& func,
                     size_t i_max,
                     vector<Res> &results,
                     Args&&... args)
  // Input: void func (size_t from, size_t to, Res& res, Args...)
  // Optimial thread_num minimizes (Time_SingleCPU/thread_num + Time_openCloseThread * (thread_num - 1)), which is sqrt(Time_SingleCPU/Time_openCloseThread)
  {
  	ASSERT (threads_max >= 1);
		results. clear ();
		results. reserve (threads_max);
  	if (threads_max == 1)
  	{
  		results. push_back (Res ());
    	func (0, i_max, results. front (), forward<Args>(args)...);
  		return;
  	}
		size_t chunk = max<size_t> (1, i_max / threads_max);
		if (chunk * threads_max < i_max)
			chunk++;
		ASSERT (chunk * threads_max >= i_max);
		Threads th (threads_max - 1);
		FFOR (size_t, tn, threads_max)
	  {
	    const size_t from = tn * chunk;
	  	if (from >= i_max)
	  		break;
	    const size_t to = from + chunk;
	    results. push_back (Res ());
	    Res& res = results. back ();
	    if (to >= i_max)
	    {
	    	func (from, i_max, res, forward<Args>(args)...);
	    	break;
	    }
		  th << thread (func, from, to, ref (res), forward<Args>(args)...);
		}
  }



//

class Verbose
{
	int verbose_old;
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
  



struct Json;
struct JsonContainer;



class Notype {};



struct Root
{
protected:
  Root () = default;
public:
  virtual ~Root () 
    {}
    // A desrtructor should be virtual to be automatically invoked by a descendant class destructor
  virtual Root* copy () const
    { NOT_IMPLEMENTED; return nullptr; }
    // Return: the same type    
  virtual void qc () const
    {}
    // Input: qc_on
  virtual void saveText (ostream& /*os*/) const 
    { NOT_IMPLEMENTED; }
    // Parsable output
  void saveFile (const string &fName) const;
    // if fName.empty() then do nothing
    // Invokes: saveText()
  string str () const
    { ostringstream oss;
      saveText (oss);
      return oss. str ();
    }
  virtual void print (ostream& os) const
    { saveText (os); }
    // Human-friendly
  virtual Json* toJson (JsonContainer* /*parent_arg*/,
                        const string& /*name_arg*/) const
    { NOT_IMPLEMENTED; return nullptr; }
	virtual bool empty () const
	  { NOT_IMPLEMENTED; return false; }
  virtual void clear ()
    { NOT_IMPLEMENTED; }
    // Postcondition: empty()
  virtual void read (istream &/*is*/)
	  { NOT_IMPLEMENTED; }
	  // Input: a line of is
};

 

inline ostream& operator<< (ostream &os,
                            const Root &r) 
  { r. saveText (os);
    return os;
  }
  


template <typename T /*Root*/> 
  struct AutoPtr : unique_ptr<T>
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
  	AutoPtr (AutoPtr<T> &t) 
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
  	AutoPtr<T>& operator= (AutoPtr<T> &t)
  	  { P::reset (t. release ());
  	  	return *this;
  	  }
  	T* copy () const
  	  { return P::get () ? static_cast<T*> (P::get () -> copy ()) : nullptr; }
  };



struct Named : Root
{
  string name;
    // !empty(), no spaces at the ends, printable ASCII characeters

  Named () = default;
  explicit Named (const string &name_arg);   
  Named* copy () const override
    { return new Named (*this); }    
  void qc () const override;
  void saveText (ostream& os) const override
    { os << name; }
	bool empty () const override
	  { return name. empty (); }
  void clear () override
    { name. clear (); }
  void read (istream &is) override
	  { is >> name; }
};


inline string named2name (const Named* n)
  { return n ? n->name : "<anonymous>"; }



typedef int (*CompareInt) (const void*, const void*);



template <typename T>
  ostream& operator<< (ostream &os,
                       const vector<T> &ts)
    { for (const auto& t : ts)
        os << t << endl;
      return os;
    }



template <typename T>
struct Vector : vector<T>
{
private:
	typedef  vector<T>  P;
public:
  bool searchSorted {false};
	

	Vector () = default;
	explicit Vector (size_t n, 
	                 const T &value = T ())
	  : P (n, value)
	  {}
  explicit Vector (initializer_list<T> init)
    : P (init)
    {}
	Vector (const vector<T> &other) 
    : P (other)
	  {}
  Vector (const Vector<T> &other) 
    : P (other)
    , searchSorted (other. searchSorted)
    {}
	Vector<T>& operator= (const Vector<T> &other)
	  { searchSorted = other. searchSorted;
	    P::reserve (other. size ());
	    P::operator= (other); 
	    return *this;
	  }
	
	
  void reserveInc (size_t inc)
    { P::reserve (P::size () + inc); }
  bool find (const T &value,
             size_t &index) const
	  // Output: index: valid if (bool)Return
	  { for (index = 0; index < P::size (); index++)
	      if (P::operator[] (index) == value)
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
  size_t countValue (const T &value) const
    { size_t n = 0;
      for (const T& t : *this)
        if (value == t)
          n++;
      return n;
    }
  Vector<T>& operator<< (const T &value)
    { P::push_back (value);
      searchSorted = false;
    	return *this;
    }
  Vector<T>& operator<< (T &&value)
    { P::push_back (move (value));
      searchSorted = false;
    	return *this;
    }
  template <typename U/*:<T>*/>
    Vector<T>& operator<< (const vector<U> &other)
      { reserveInc (other. size ());
        P::insert (P::end (), other. begin (), other. end ());
        searchSorted = false;
      	return *this;
      }
  template <typename U/*:<T>*/>
    Vector<T>& operator<< (const list<U> &other)
      { reserveInc (other. size ());
        P::insert (P::end (), other. begin (), other. end ());
        searchSorted = false;
      	return *this;
      }
  void setAll (const T &value)
    { for (T &t : *this)
    	  t = value;
    	searchSorted = false;
    }
  void eraseAt (size_t index)
    { eraseMany (index, index + 1); }
  void eraseMany (size_t from,
                  size_t to)
    { P::erase ( P::begin () + (ptrdiff_t) from
    	         , P::begin () + (ptrdiff_t) to
    	         ); 
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
    	  swap (P::operator[] (i), P::operator[] (j));
    	}
    	searchSorted = false;
    }
  void randomOrder ()
		{ Rand rand (seed_global);
			for (T &t : *this)
	      swap (t, P::operator[] ((size_t) rand. get ((ulong) P::size ())));
    	searchSorted = false;
		}
  T pop (size_t n = 1)
    { T t = T ();
      while (n)
      { t = P::operator[] (P::size () - 1);
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
            (*this) [j] = (*this) [i];
          if (cond (j))
            toDelete++;
        }
        while (toDelete)
        { P::pop_back ();
          toDelete--;
        }
      }
  template <typename Condition /*on index*/>
    void filterValue (const Condition cond)
      { size_t toDelete = 0;
        for (size_t i = 0, end_ = P::size (); i < end_; i++)
        { const size_t j = i - toDelete;
          if (j != i)
            (*this) [j] = (*this) [i];
          if (cond ((*this) [j]))
            toDelete++;
        }
        while (toDelete)
        { P::pop_back ();
          toDelete--;
        }
      }

  void sort ()
    { if (searchSorted)
        return;
      Common_sp::sort (*this); 
      searchSorted = true;
    }
  template <typename StrictlyLess>
    void sort (const StrictlyLess &strictlyLess)
      { Common_sp::sort (*this, strictlyLess); 
        searchSorted = false;
      }    
  void sortBubble ()  
    // Fast if *this is almost sort()'ed
    { if (searchSorted)
        return;
      FOR_START (size_t, i, 1, P::size ())
		    FOR_REV (size_t, j, i)
		      if (P::operator[] (j + 1) > P::operator[] (j))
        	  swap (P::operator[] (j), P::operator[] (j + 1));
		      else
		      	break;
    	searchSorted = true;
    }
  void checkSorted () const
    { if (! searchSorted)
    	  throw logic_error ("Vector is not sorted for search");
    }

  size_t binSearch (const T &value,
                    bool exact = true) const
    // Return: if exact then NO_INDEX or vec[Return] = value else min {i : vec[i] >= value}
    { if (P::empty ())
    	  return NO_INDEX;
    	checkSorted ();
    	size_t lo = 0;  // vec.at(lo) <= value
    	size_t hi = P::size () - 1;  
    	// lo <= hi
    	if (value < P::operator[] (lo))
    	  return exact ? NO_INDEX : lo;
    	if (P::operator[] (hi) < value)
    	  return NO_INDEX;
    	// at(lo) <= value <= at(hi)
    	for (;;)
    	{
	    	const size_t m = (lo + hi) / 2;
	    	if (   P::operator[] (m) == value
	    		  || P::operator[] (m) <  value
	    		 )
	    		if (lo == m)  // hi in {lo, lo + 1}
	    			break;
	    		else
	    		  lo = m;
	      else
	      	hi = m;
	    }
	    if (P::operator[] (lo) == value)
	    	return lo;
	    if (! exact || P::operator[] (hi) == value)
	    	return hi;
	    return NO_INDEX;
    }
  template <typename U /* : T */>
    bool containsFast (const U &value) const
      { return binSearch (value) != NO_INDEX; }
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
    bool intersectsFast (const Vector<U> &other) const
      { for (const U& u : other)
          if (containsFast (u))
            return true;
        return false;
      }
  template <typename U>
    bool intersectsFast_merge (const Vector<U> &other) const
      { checkSorted ();
      	other. checkSorted ();
      	size_t i = 0;
      	const size_t otherSize = other. size ();
        for (const T& t : *this)
        { while (i < otherSize && other [i] < t)
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
      { for (const T& t : *this)
          if (other. find (t) != other. end ())
            return true;
        return false;
      }
  template <typename U /* : T */>
    void setMinus (const Vector<U> &other)
      { filterIndex ([&] (size_t i) { return other. containsFast ((*this) [i]); }); }
      
  size_t findDuplicate () const
    { if (P::size () <= 1)
        return NO_INDEX;
      FOR_START (size_t, i, 1, P::size ())
        if ((*this) [i] == (*this) [i - 1])
          return i;
      return NO_INDEX;
    }
  bool isUniq () const
    { return findDuplicate () == NO_INDEX; }
  template <typename Equal /*bool equal (const T &a, const T &b)*/>
	  void uniq (const Equal &equal)
	    { if (P::size () <= 1)
	        return;
	      size_t j = 1;  
	      FOR_START (size_t, i, 1, P::size ())
	        if (! equal ((*this) [i], (*this) [i - 1]))
	        { if (j != i)
	            (*this) [j] = (*this) [i];
	          j++;
	        }
	      P::resize (j);
	    }
	void uniq ()
	  { uniq ([] (const T& a, const T& b) { return a == b; }); }
  size_t getIntersectSize (const Vector<T> &other) const
    // Input: *this, vec: unique
    { if (other. empty ())
        return 0;
      checkSorted ();
      other. checkSorted ();      
      size_t n = 0;
      size_t j = 0;
      for (const T& x : *this)
      { while (other [j] < x)
        { j++;
          if (j == other. size ())
            return n;
        }
        if (other [j] == x)
          n++;
      }
      return n;
    }

  bool operator< (const Vector<T> &other) const
    // Lexicographic comparison
    { FFOR (size_t, i, std::min (P::size (), other. size ()))
    	{ if (P::operator[] (i) < other [i]) return true;
    		if (other [i] < P::operator[] (i)) return false;
      }
      return P::size () < other. size ();
    }
};



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
	VectorPtr (const VectorPtr<T> & /*other*/) = default;
  explicit VectorPtr (initializer_list<const T*> init)
    : P (init)
    {}
	template <typename U>
  	VectorPtr (const vector<const U*> &other)
  	  : P ()
  	  { P::reserve (other. size ());
  	    insertAll (*this, other);
  	  }	  
	template <typename U>
  	VectorPtr (const list<const U*> &other)
  	  : P ()
  	  { P::reserve (other. size ());
  	    insertAll (*this, other);
  	  }	  


	void deleteData ()
	  {	for (const T* t : *this)
			  delete t;
			P::clear ();  
	  }
  void erasePtr (size_t index)
    { delete P::operator[] (index);
      P::eraseAt (index);
    }
  void sortBubblePtr ()
    { FOR_START (size_t, i, 1, P::size ())
		    FOR_REV (size_t, j, i)
		      if (* P::operator[] (j + 1) > * P::operator[] (j))
        	  std::swap (P::operator[] (j), P::operator[] (j + 1));
		      else
		      	break;
    }
};



template <typename T /* : Root */>
struct VectorOwn : VectorPtr<T>
{
private:
	typedef  VectorPtr<T>  P;
public:

  VectorOwn () = default;
	VectorOwn (const VectorOwn<T> &x)
	  : P ()
	  { *this = x; }
	VectorOwn<T>& operator= (const VectorOwn<T> &x)
	  { P::deleteData ();
	  	P::reserve (x. size ());
	  	for (const T* t : x)
	  	  P::push_back (static_cast <const T*> (t->copy ()));
	  	P::searchSorted = x. searchSorted;
	  	return *this;
	  }
	VectorOwn (const VectorPtr<T> &x)
	  : P ()
	  { P::operator= (x); }
 ~VectorOwn ()
    { P::deleteData (); }
};



struct StringVector : Vector<string>
{
private:
	typedef  Vector<string>  P;
public:
	

  StringVector () = default;
  StringVector (initializer_list<string> init)
    : P (init)
    {}
  StringVector (const string &fName,
                size_t reserve_size);
  explicit StringVector (const string &s, 
                         char c = ' ');


  string toString (const string& sep) const
    { string res;
  	  for (const string& s : *this)
  	  { if (! res. empty ())
  	      res += sep;
  	    res += s;
  	  }
  	  return res;
  	}
};



struct DisjointCluster
// Cormen, Leiserson, Rivest, Introduction to Algorithms, p. 449
{
protected:
  DisjointCluster* parentDC;
    // !nullptr
    // Tree
    // = this <=> root
  size_t rankDC;
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



template <typename T>
struct Heap : Root, Nocopy
// Priority queue
// Heap property: comp(&arr[parent(index)],&arr[index]) >= 0
// More operations than in std::priority_queue
{
private:
  Vector<T> arr;
    // Elements are not owned by arr
  const CompareInt comp;
    // !nullptr
  typedef void (*SetHeapIndex) (T &item, size_t index);
    // Example: item.heapIndex = index
  const SetHeapIndex setHeapIndex;
    // Needed to invoke increaseKey()
public:


  struct Error : runtime_error
    { explicit Error (const string &str) : runtime_error (("Heap: " + str)) {} };


  explicit Heap (const CompareInt &comp_arg,
					       const SetHeapIndex &setHeapIndex_arg = nullptr,
					       size_t toReserve = 0)
    : comp (comp_arg)
    , setHeapIndex (setHeapIndex_arg)
    { arr. reserve (toReserve); }


  bool empty () const final
    { return arr. empty (); }
  Heap& operator<< (T item)
    { arr << item;
      increaseKey (arr. size () - 1);
      return *this;
    }
  T increaseKey (size_t index)
    { T item = arr [index];
      size_t p;
      while (index && comp (& arr [p = parent (index)], & item) < 0)
      { assign (arr [p], index);
        index = p;
      }
      assign (item, index);
      return item;
    }
  T decreaseKey (size_t index)
    { T item = arr [index];
      heapify (index, arr. size ());
      return item;
    }
  T getMaximum () const
    { if (arr. empty ()) 
    	  throw Error ("getMaximum");
      return arr [0];
    }
  void deleteMaximum ()
    // Time: O(1) amortized
    { if (arr. empty ()) 
    	  throw Error ("deleteMaximum");
      T item = arr. back ();
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
  size_t parent (size_t index) const
    { if (! index) throw Error ("parent");
      return (index + 1) / 2 - 1;
    }
  size_t left (size_t index) const
    { return 2 * (index + 1) - 1; }
  size_t right (size_t index) const
    { return left (index) + 1; }
  void assign (T item,
               size_t index)
    { arr [index] = item;
      if (setHeapIndex)
        setHeapIndex (item, index);
    }
  void swap (size_t i,
             size_t j)
    { T item = arr [i];
      assign (arr [j], i);
      assign (item, j);
    }
  void heapify (size_t index,
                size_t maxIndex)
    // Requires: Heap property holds for all index1 < maxIndex except parent(index1) == index
    { if (maxIndex > arr. size ()) throw Error ("heapify: maxIndex");
      if (index >= maxIndex)       throw Error ("heapify: index");
      for (;;)
      { size_t extr = index;
        const size_t l = left (index);
        if (   l < maxIndex
            && comp (& arr [extr], & arr [l]) < 0)
          extr = l;
        const size_t r = right (index);
        if (   r < maxIndex
            && comp (& arr [extr], & arr [r]) < 0)
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
    { Heap <string> heap (strComp);
      heap << "Moscow" << "San Diego" << "Los Angeles" << "Paris";
      while (! heap. empty ())  
      { cout << heap. getMaximum () << endl;
        heap. deleteMaximum ();
      }
    }
private:
  static int strComp (const void* s1,
                      const void* s2)
    { const string& s1_ = * static_cast <const string*> (s1);
      const string& s2_ = * static_cast <const string*> (s2);
      if (s1_ < s2_) return -1;
      if (s1_ > s2_) return  1;
      return  0;
    }
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
	Set (const Set<T> &other)
	  : P ()
	  { *this = other; }
	Set<T>& operator= (const Set<T> &other)
	  { universal = other. universal;
	  	return static_cast <Set<T>&> (P::operator= (other)); 
	  }
	template <typename U, typename V>
	  Set (const map<U,V> &other)
	    : universal (false)
	    { for (const auto& it : other)
	        P::insert (it. first);
	    }
	template <typename U, typename V>
	  Set (const unordered_map<U,V> &other)
	    : universal (false)
	    { for (const auto& it : other)
	        P::insert (it. first);
	    }
	template <typename U>
	  Set (const vector<U> &other)
	    : universal (false)
	    { for (const U& u : other)
	        P::insert (u);
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
    { for (const T& t : *this)
        if (other. contains (t))
          return true;
      return false;
    }
public:

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
      for (Iter <Set<T>> iter (*this); iter. next (); )
				if (! other. contains (*iter))
					iter. erase ();
		}
	size_t intersectSize (const Set<T> &other) const
	  // Return: universal <=> SIZE_MAX
		{ if (other. universal)
			  return universal ? SIZE_MAX : P::size ();
			if (universal)
				return other. size ();
		  size_t n = 0;
		  for (const T& t : *this)
				if (other. contains (t))
					n++;
			return n;
		}
  size_t setMinus (const Set<T> &other)
    { ASSERT (! universal);
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
    	IMPLY (from, ! from->universal);
    	IMPLY (to,   ! to  ->universal);
    	if (from)
    	  { EXEC_ASSERT (from->erase (el) == 1); }
    	if (to)
    	  { EXEC_ASSERT (to->insert (el). second); }
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



struct Progress : Nocopy
{
private:
	static size_t beingUsed;
public:

	size_t n_max {0};
	  // 0 <=> unknown
	bool active;
	size_t n {0};
	string step;
	size_t displayPeriod {0};
	

	explicit Progress (size_t n_max_arg = 0,
	                   size_t displayPeriod_arg = 1)
	  : n_max (n_max_arg)
	  , active (enabled () && displayPeriod_arg && (! n_max_arg || displayPeriod_arg <= n_max_arg))
	  , displayPeriod (displayPeriod_arg)
	  { if (active) 
	  	  beingUsed++; 
	  	if (active)
	  	  report ();
	  }
 ~Progress () 
    { if (active)
    	{ if (! uncaught_exception ())
    	  { report ();
    	    cerr << endl;
    	  }
    	  beingUsed--;
    	}
    }
    

  void operator() (const string& step_arg = string ())
    { n++;
    	step = step_arg;
    	if (   active 
    		  && n % displayPeriod == 0
    		 )
    	  report ();
    }
private:
	void report () const;
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
	  { return ! beingUsed && verbose (1); }
};




// Input

struct Input : Root, Nocopy
{
protected:
	AutoPtr <char> buf;
	ifstream ifs;
  istream* is {nullptr};
    // !nullptr
public:
	bool eof {false};
	uint lineNum {0};
	  // # lines read
protected:
	Progress prog;
public:


protected:	
  Input (const string &fName,
         size_t bufSize,
         uint displayPeriod);
  Input (istream &is_arg,
	       uint displayPeriod);
public:


  void reset ();
    // Update: ifs
};
	


struct LineInput : Input
{
	string line;
	  // Current line
  string commentStart;

	
	explicit LineInput (const string &fName,
          	          size_t bufSize = 100 * 1024,
          	          uint displayPeriod = 0)
    : Input (fName, bufSize, displayPeriod)
    {}
  explicit LineInput (istream &is_arg,
	                    uint displayPeriod = 0)
    : Input (is_arg, displayPeriod)
    {}


	bool nextLine ();
  	// Output: eof, line
  	// Update: lineNum
    // Invokes: trimTrailing()
	bool expectPrefix (const string &prefix,
	                   bool eofAllowed)
		{ if (nextLine () && trimPrefix (line, prefix))
		  	return true;  
			if (eof && eofAllowed)
				return false;
		  throw runtime_error ("No " + strQuote (prefix));
		}
	string getString ()
	  { string s; 
	  	while (nextLine ())
	  	{ if (! s. empty ())
	  			s += "\n";
	  	  s += line;
	  	}
	  	return s;
	  }
	StringVector getVector ()
	  { StringVector vec;
	    while (nextLine ())
	      if (! line. empty ())
	  	    vec << line;
	  	return vec;
	  }
};
	


struct ObjectInput : Input
{
	explicit ObjectInput (const string &fName,
          	            size_t bufSize = 100 * 1024,
          	            uint displayPeriod = 0)
    : Input (fName, bufSize, displayPeriod)
    {}
  explicit ObjectInput (istream &is_arg,
	                      uint displayPeriod = 0)
    : Input (is_arg, displayPeriod)
    {}


	bool next (Root &row);
	  // Output: row
  	// Update: lineNum
};
	


struct CharInput : Input
{
  uint charNum {(uint) -1};
    // In the current line
  bool eol {false};
    // eof => eol
private:
  bool ungot {false};
public:

	
	explicit CharInput (const string &fName,
              	      size_t bufSize = 100 * 1024,
              	      uint displayPeriod = 0)
    : Input (fName, bufSize, displayPeriod)
    {}
  explicit CharInput (istream &is_arg,
	                    uint displayPeriod = 0)
    : Input (is_arg, displayPeriod)
    {}


	char get ();
	  // Output: eof
	  // Update: lineNum, charNum
	void unget ();
	  // Requires: To be followed by get()
  string getLine ();
    // Postcondition: eol
	  

  struct Error : runtime_error
    { explicit Error (const CharInput &in,
		                  const string &expected = string ()) 
			  : runtime_error ("Error at line " + toString (in. lineNum + 1) 
		                     + ", pos. " + toString (in. charNum + 1)
		                     + (expected. empty () ? string () : (": " + expected + " is expected"))
		                    )
	      {}
	  };
};
	


struct PairFile : Root
{
private:
	LineInput f;
public:
	string name1;
	string name2;
	// name1 < name2
	
	explicit PairFile (const string &fName)
	  : f (fName, 100 * 1024, 1000)  // PAR
	  {}
	  
	bool next ()
	  { if (! f. nextLine ())
	  	  return false;
	  	static Istringstream iss;
      iss. reset (f. line);
      iss >> name1 >> name2;
      if (name2. empty ())
      	throw runtime_error ("Bad request: " + strQuote (name1) + " - " + strQuote (name2));
      if (name1 == name2)
      	throw runtime_error ("Same name: " + name1);
      if (name1 > name2)
      	swap (name1, name2);
      return true;
	  }
};



struct Token : Root
{
	static constexpr char quote = '\"';
	static constexpr uint noNum = (uint) -1;
	enum Type { eText       // In quote's within one line
	          , eNumber     // All characters: isDigit()
	          , eName       // All characters: isLetter()
	                        // Start: !isDigit()
	          , eDelimiter  // Does not include ' '
	          };
	Type type {eText};
	string name;
	  // All characters: printable()
	  // eText => embracing quote's are removed
	  // May be: !goodName(name)
	uint num {noNum};
	  // Valid if type = eNumber
	uint charNum {0};	
	  // First CharInput::charNum of name

	  
	Token ()
	  { clear (); }
	Token (const string& name_arg,
	       Type type_arg = eName)
	  : type (type_arg)
	  , name (name_arg)
	  {}
	Token (const uint num_arg)
	  : type (eNumber)
	  , num (num_arg)
	  {}
	Token (char delimiter_arg)
	  : type (eDelimiter)
	  , name (1, delimiter_arg)
	  {}
	explicit Token (CharInput &in)
	  { readInput (in); }
	Token (CharInput &in,
	       Type expected)
    { readInput (in);
    	if (empty ())
 			  throw CharInput::Error (in, "No token"); 
    	if (type != expected)
 			  throw CharInput::Error (in, type2str (type)); 
    }
	Token (CharInput &in,
	       uint expected)
    { readInput (in);
    	if (! isNumber (expected))
 			  throw CharInput::Error (in, type2str (eNumber) + " " + toString (expected)); 
    }
	Token (CharInput &in,
	       const string &expected)
    { readInput (in);
    	if (! isName (expected))
 			  throw CharInput::Error (in, type2str (eName) + " " + expected); 
    }
	Token (CharInput &in,
	       char expected)
    { readInput (in);
    	if (! isDelimiter (expected))
 			  throw CharInput::Error (in, type2str (eDelimiter) + " " + expected); 
    }
private:
	void readInput (CharInput &in);
public:
	void qc () const override;
	void saveText (ostream &os) const override
	  { if (! empty ())
  	    switch (type)
  	  	{ case eText:      os << ' ' << quote << name << quote; break;
  	  		case eNumber:    os << ' ' << num;                    break;
  	  		case eName:      os << ' '          << name;          break;
  	  		case eDelimiter: os                 << name;          break;
  	  	}
	  }
	bool empty () const override
	  { return type == eDelimiter && name. empty (); }
	void clear () override
	  { type = eDelimiter;
	    name. clear ();
	  	num = noNum;
	  	charNum = 0;
	  }


	static string type2str (Type type) 
	  { switch (type)
	  	{ case eText:      return "Text";
	  		case eNumber:    return "Number";
	  		case eName:      return "Name";
	  		case eDelimiter: return "Delimiter";
	  	}
	  	return "?";
	  }
	bool isNumber (uint n) const
	  { return ! empty () && type == eNumber && num == n; }
	bool isName (const string &s) const
	  { return ! empty () && type == eName && name == s; }
	bool isDelimiter (char c) const
	  { return ! empty () && type == eDelimiter && name [0] == c; }
};




struct OFStream : ofstream
{
	OFStream () = default;
	OFStream (const string &dirName,
	          const string &fileName,
	          const string &extension)
	  { open (dirName, fileName, extension); }
	explicit OFStream (const string &pathName)
	  { open ("", pathName, ""); }


	void open (const string &dirName,
	           const string &fileName,
	           const string &extension);
	  // Input: !fileName.empty()
};



struct Stderr : Singleton<Stderr>
{
  bool quiet {false};
  
  explicit Stderr (bool quiet_arg)
    : quiet (quiet_arg)
    {}
    
  template <typename T>
    Stderr& operator<< (const T& t) 
      { if (quiet)
          return *this;
        cerr << t;
        return *this;
      }
};



struct Csv : Root
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

  
  
void csvLine2vec (const string &line,
                  StringVector &words);
  // Output: words
  // Invokes: Csv




struct TabDel
// Usage: {<<field;}* str();
{
private:
  ostringstream tabDel;
  ONumber on;
public:
  
  TabDel (streamsize precision = 6,
	        bool scientific = false)
	  : on (tabDel, precision, scientific)
	  {}
    
  template <typename T>
    TabDel& operator<< (const T &field)
      { if (tabDel. tellp ())  
          tabDel << '\t'; 
        tabDel << field; 
        return *this; 
      }    
  string str () const
    { return tabDel. str (); }
};




// Json

struct JsonNull;
struct JsonInt;
struct JsonDouble;
struct JsonString;
struct JsonBoolean;
struct JsonArray;
struct JsonMap;
  

struct Json : Root, Nocopy  // Heaponly
{
protected:
  Json (JsonContainer* parent,
        const string& name);
  Json () = default;
public:  
  void print (ostream& os) const override = 0;
  
  virtual const JsonNull* asJsonNull () const
    { return nullptr; }  
  virtual const JsonInt* asJsonInt () const
    { return nullptr; }  
  virtual const JsonDouble* asJsonDouble () const
    { return nullptr; }  
  virtual const JsonString* asJsonString () const
    { return nullptr; }  
  virtual const JsonBoolean* asJsonBoolean () const
    { return nullptr; }  
  virtual const JsonArray* asJsonArray () const
    { return nullptr; }  
  virtual const JsonMap* asJsonMap () const
    { return nullptr; }  

protected:
  static string toStr (const string& s)
    { return "'" + to_c (s) + "'"; }    
  static Token readToken (istream &is);
  static void parse (istream& is,
                     const Token& firstToken,
                     JsonContainer* parent,
                     const string& name);
public:
    
  int getInt () const;
    // Requires: JsonInt
  double getDouble () const;
    // Requires: JsonDouble
  string getString () const;
    // Requires: JsonString
  const Json* at (const string& name_arg) const;
    // Requires: JsonMap
  const Json* at (size_t index) const;
    // Requires: JsonArray
  size_t getSize () const;
    // Requires: JsonArray
};


struct JsonNull : Json
{
  explicit JsonNull (JsonContainer* parent,
                     const string& name = noString)
    : Json (parent, name)
    {}    
  void print (ostream& os) const final
    { os << "null"; }

  const JsonNull* asJsonNull () const final
    { return this; }  
};


struct JsonInt : Json
{
  int n;
  
  JsonInt (int n_arg,
           JsonContainer* parent,
           const string& name = noString)
    : Json (parent, name)
    , n (n_arg)
    {}
  void print (ostream& os) const final
    { os << n; }

  const JsonInt* asJsonInt () const final
    { return this; }  
};


struct JsonDouble : Json
{
  double n;
  streamsize decimals;

  JsonDouble (double n_arg,
              streamsize decimals_arg,
              JsonContainer* parent,
              const string& name = noString)
    : Json (parent, name)
    , n (n_arg)
    , decimals (decimals_arg == numeric_limits<streamsize>::max() ? double2decimals (n_arg) : decimals_arg)
    {}
    // decimals_arg = -1: default
  void print (ostream& os) const final
    { const ONumber on (os, (streamsize) decimals, false);
    	if (n == n)
        os << n; 
      else
        os << "null";  // NaN
    }      

  const JsonDouble* asJsonDouble () const final
    { return this; }  
};


struct JsonString : Json
{
  string s;

  JsonString (const string& s_arg,
              JsonContainer* parent,
              const string& name = noString)
    : Json (parent, name)
    , s (s_arg)
    {}
  void print (ostream& os) const final
    { os << toStr (s); }

  const JsonString* asJsonString () const final
    { return this; }  
};


struct JsonBoolean : Json
{
  bool b;

  JsonBoolean (bool b_arg,
               JsonContainer* parent,
               const string& name = noString)
    : Json (parent, name)
    , b (b_arg)
    {}
  void print (ostream& os) const final
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


struct JsonArray : JsonContainer
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
  JsonArray (istream& is,
             JsonContainer* parent,
             const string& name);
public:
  void print (ostream& os) const final;

  const JsonArray* asJsonArray () const final
    { return this; }
};


struct JsonMap : JsonContainer
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
  JsonMap (istream& is,
           JsonContainer* parent,
           const string& name)
    : JsonContainer (parent, name)
    { parse (is); }
  void parse (istream& is);
public:
 ~JsonMap ();
  void print (ostream& os) const final;

  const JsonMap* asJsonMap () const final
    { return this; }
};



extern JsonMap* jRoot;



//

struct Offset 
// Linked to the ostream os
// Not thread-safe
{
private:
	static size_t size;
public:
	static constexpr size_t delta = 2;

	Offset ()
  	{ size += delta; }
 ~Offset () 
  	{ size -= delta; }

  static void newLn (ostream &os) 
    { os << endl << string (size, ' '); }
};



void exec (const string &cmd,
           const string &logFName = string());




///////////////////////////////////////////////////////////////////////////

struct ItemGenerator
{
  Progress prog;
  
  ItemGenerator (size_t progress_n_max,
	               size_t progress_displayPeriod)
	  : prog (progress_n_max, progress_displayPeriod)
	  {}
  virtual ~ItemGenerator ()
    {}
  
  virtual bool next (string &item) = 0;
    // Return: false <=> end of items
    // Output: item; may be empty()
    // Invokes: prog()
};



struct FileItemGenerator : ItemGenerator, Nocopy
{
  const bool isDir;
private:
  string fName;
  ifstream f;
public:
  
  FileItemGenerator (size_t progress_displayPeriod,
                     bool isDir_arg,
                     const string& fName_arg);
 ~FileItemGenerator ()
    { if (isDir)
	      remove (fName. c_str ());
	  }
  
  bool next (string &item) final;
};

  

struct NumberItemGenerator : ItemGenerator
{
private:
  const size_t n;
  size_t i {0};
public:
  
  NumberItemGenerator (size_t progress_displayPeriod,
                       const string& name)
    : ItemGenerator (str2<size_t> (name), progress_displayPeriod)
    , n (prog. n_max)
    {}
  
  bool next (string &item) final
    { if (i == n)
        return false;
      i++;
      item = toString (i);
      prog (item);
      return true;
    }
};

  


///////////////////////////////////////////////////////////////////////////

struct Application : Singleton<Application>, Root
// Usage: int main (argc, argv) { ThisApplication /*:Application*/ app; return app. run (argc, argv); }
{  
  const string description;
  const bool needsArg;
  const bool gnu;
  string version {"1"};
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
  struct Positional : Arg
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
  struct Key : Arg
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
  List<Positional> positionals;
  List<Key> keys;
  map<string/*Arg::name*/,const Arg*> name2arg;
  map<char/*Arg::name[0]*/,const Key*> char2arg;
    // Valid if gnu
public:
  

protected:
  explicit Application (const string &description_arg,
                        bool needsArg_arg = true,
                        bool gnu_arg = false)               
    : description (description_arg)
    , needsArg (needsArg_arg)
    , gnu (gnu_arg)
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
               const string &defaultValue = string (),
               char acronym = '\0',
               const string &var = string ());
    // [-<name> <defaultValue>]
    // gnu: [--<name> <var>]
  void addFlag (const string &name,
                const string &argDescription,
                char acronym = '\0');
    // [-<name>]
  void addPositional (const string &name,
                      const string &argDescription);
  void setRequiredGroup (const string &keyName,
                         const string &requiredGroup);
private:
	void addDefaultArgs ()
	  { if (gnu)
    	{ addKey ("threads", "Max. number of threads", "1", '\0', "THREADS");
    	  addFlag ("debug", "Integrity checks");
      }
    	else
    	{ addFlag ("qc", "Integrity checks (quality control)");
	      addKey ("verbose", "Level of verbosity", "0");
	      addFlag ("noprogress", "Turn off progress printout");
	      addFlag ("profile", "Use chronometers to profile");
	      addKey ("seed", "Positive integer seed for random number generator", "1");
	      addKey ("threads", "Max. number of threads", "1");
	      addKey ("json", "Output file in Json format");
	      addKey ("log", "Error log file, appended");
	      addFlag ("sigpipe", "Exit normally on SIGPIPE");
	    }
	  }
	void qc () const final;
	Key* getKey (const string &keyName) const;
	  // Return: !nullptr
	void setPositional (List<Positional>::iterator &posIt,
	                    const string &value);
public:
  virtual ~Application ()
    {}


protected:
  string getArg (const string &name) const;
    // Input: keys, where Key::flag = false, and positionals
  uint arg2uint (const string &name) const
    { uint n = 0;
    	try { n = str2<uint> (getArg (name)); }
    	  catch (...) { throw runtime_error ("Cannot convert -" + name + " to non-negative number"); }
    	return n;
    }
  double arg2double (const string &name) const
    { double d = numeric_limits<double>::quiet_NaN ();
    	try { d = str2<double> (getArg (name)); }
    	  catch (...) { throw runtime_error ("Cannot convert -" + name + " to number"); }
    	return d;
    }
  bool getFlag (const string &name) const;
    // Input: keys, where Key::flag = true
  string key2shortHelp (const string &name) const;
  string getProgramDirName () const
    { return getDirName (programArgs. front ()); }
protected:
  virtual void initEnvironment ()
    {}
  string getInstruction () const;
  string getHelp () const;
public:
  int run (int argc, 
           const char* argv []);
    // Invokes: body()
private:
  virtual void body () const = 0;
    // Invokes: initEnvironment()
};



struct ShellApplication : Application
// Requires: SHELL=bash
{
  // Environment
  const bool useTmp;
  string tmp;
  string execDir;
  mutable map<string,string> prog2dir;
  

  ShellApplication (const string &description_arg,
                    bool needsArg_arg,
                    bool gnu_arg,
                    bool useTmp_arg)
    : Application (description_arg, needsArg_arg, gnu_arg)
    , useTmp (useTmp_arg)
    {}
 ~ShellApplication ();


protected:
  void initEnvironment () override;
private:
  void body () const final;
  virtual void shellBody () const = 0;
protected:
  static string shellQuote (string s)
    { replaceStr (s, "\'", "\'\"\'\"\'");
    	return "\'" + s + "\'";
    }
  static bool emptyArg (const string &s)
    {	return s. empty () || s == "\'\'"; }
  string which (const string &progName) const;
    // Return: isRight(,"/") or empty()
  void findProg (const string &progName) const;
    // Output: prog2dir
  string fullProg (const string &progName) const;
    // Return: directory + progName + ' '
    // Requires: After findProg(progName)
};



}



#endif

