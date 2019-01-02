// common.cpp

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


#undef NDEBUG
#include "common.inc"

#include "common.hpp"

#include <sstream>
#include <cstring>
#ifndef _MSC_VER
  #include <dirent.h>
#endif
#include <signal.h>  



namespace Common_sp
{
 


vector<string> programArgs;
string programName;



string getCommandLine ()
{
  string commandLine;
  for (const string& s : programArgs)
  {
    const bool bad =    s. empty () 
                     || contains (s, ' ')
                     || contains (s, '|')
                     || contains (s, ';')
                     || contains (s, '#')
                     || contains (s, '*')
                     || contains (s, '?')
                     || contains (s, '$')
                     || contains (s, '(')
                     || contains (s, ')')
                     || contains (s, '<')
                     || contains (s, '>')
                     || contains (s, '~')
                      || contains (s, '\'')
                     || contains (s, '"')
                     || contains (s, '\\');
    if (! commandLine. empty ())
      commandLine += ' ';
    if (bad)
      commandLine += strQuote (s, '\'');
    else
      commandLine += s;
  }
  return commandLine;
}



ostream* logPtr = nullptr;

bool qc_on = false;
ulong seed_global = 1;
bool sigpipe = false;


// thread
size_t threads_max = 1;

inline thread::id get_thread_id ()
  { return this_thread::get_id (); }

thread::id main_thread_id = get_thread_id ();
  
bool isMainThread ()
  { return get_thread_id () == main_thread_id; }


bool Chronometer::enabled = false;



namespace 
{

void segmFaultHandler (int /*sig_num*/)
{
  signal (SIGSEGV, SIG_DFL); 
  errorExit ("Segmentation fault", true); 
}


void sigpipe_handler (int /*sig_num*/) 
{
  signal (SIGPIPE, SIG_DFL);
  if (sigpipe)
    exit (0);
  exit (1);
}

}



bool initCommon ()
{
  MODULE_INIT

  signal (SIGSEGV, segmFaultHandler);  
  signal (SIGPIPE, sigpipe_handler);
        
#ifdef _MSC_VER
  #pragma warning (disable : 4127)
#endif
  ASSERT (numeric_limits<char>::is_signed);
  ASSERT (sizeof (long) >= 4);
#ifdef _MSC_VER
  #pragma warning (default : 4127)
#endif

  ASSERT (SIZE_MAX == std::numeric_limits<size_t>::max ());

  return true;
}


namespace
{
  const bool init_ = initCommon ();
}



void errorExit (const char* msg,
                bool segmFault)
// alloc() may not work
{ 
	ostream* os = logPtr ? logPtr : & cout; 
	
	// time ??
#ifndef _MSC_VER
	const char* hostname = getenv ("HOSTNAME");
	const char* pwd = getenv ("PWD");
#endif
  const string errorS ("*** ERROR ***");
  if (isLeft (msg, errorS))  // msg already is the result of errorExit()
    *os << endl << msg << endl;
  else
  	*os << endl
        << errorS << endl
        << msg << endl << endl
      #ifndef _MSC_VER
  	    << "HOSTNAME: " << (hostname ? hostname : "?") << endl
  	    << "PWD: " << (pwd ? pwd : "?") << endl
      #endif
  	    << "Progam name:  " << programName << endl
  	    << "Command line: " << getCommandLine () << endl;
  //system (("env >> " + logFName). c_str ());

  os->flush ();

  if (segmFault)
    abort ();
  exit (1);
}



namespace
{
	
uint powInt_ (uint a,
              uint b)
// Input: a: !0, != 1
{
	if (! b)
		return 1;
	if (b == 1)
		return a;
	const uint half = b / 2;
	const uint res = powInt_ (a, half);
	return res * res * (divisible (b, 2) ? 1 : a);
}
	
}


uint powInt (uint a,
             uint b)
{
	if (a)
		if (a == 1)
			return 1;
		else
			return powInt_ (a, b);
	else
		if (b)
			return 0;
		else
			throw runtime_error ("powInt: 0^0");
}




// string

const string noString;



bool isRight (const string &s,
              const string &right)
{
  if (s. size () < right. size ())
    return false;
  return s. substr (s. size () - right. size ()) == right;
}



bool trimPrefix (string &s,
                 const string &prefix)
{
  if (! isLeft (s, prefix))
    return false;
  s. erase (0, prefix. size ());
  return true;
}



bool trimSuffix (string &s,
                 const string &suffix)
{
  if (! isRight (s, suffix))
    return false;
  s. erase (s. size () - suffix. size ());
  return true;
}



void trimSuffixNonAlphaNum (string &s)
{
  for (;;)
  {
    trimTrailing (s);
    if (s. empty ())
      break;
    const size_t pos = s. size () - 1;
    if (   isAlpha (s [pos])
        || isDigit (s [pos])
       )
      break;
    s. erase (pos);
  }
}



bool trimTailAt (string &s,
                 const string &tailStart)
{
  const size_t pos = s. find (tailStart);
  const bool trimmed = pos != string::npos;
  if (trimmed)
    s. erase (pos);
  return trimmed;
}



bool goodName (const string &name)
{
  if (name. empty ())
    return false;
  if (name [0] == ' ')
    return false;
  if (*(name. end () - 1) == ' ')
    return false;

  for (const char c : name)
    if (! printable (c))
      return false;
      
  return true;
}



bool strBlank (const string &s)
{
  for (const char c : s)
    if (! isSpace (c))
      return false;
  return true;
}



void strUpper (string &s)
{
  for (char& c : s)
    c = toUpper (c);
}



void strLower (string &s)
{
  for (char& c : s)
    c = toLower (c);
}



bool isUpper (const string &s)
{
  for (const char c : s)
    if (! isUpper (c))
    	return false;
  return true;
}



bool isLower (const string &s)
{
  for (const char c : s)
    if (! isLower (c))
    	return false;
  return true;
}



string::const_iterator stringInSet (const string &s,
                    	           	  const string &charSet)
{
  CONST_ITER (string, it, s)
    if (! charInSet (*it, charSet))
      return it;

  return s. end ();
}



size_t strCountSet (const string &s,
 		                const string &charSet)
{
  size_t n = 0;
  for (const char c : s)
	  if (charInSet (c, charSet))		
	    n++;
  return n;
}



void strDeleteSet (string &s,
		               const string &charSet)
{
  FOR_REV (size_t, i, s. size ())
    if (charInSet (s [i], charSet))
    	s. erase (i, 1);
}



void trimLeading (string &s)
{
  size_t i = 0;
  for (; i < s. size () && isSpace (s [i]); i++)
    ;
  s. erase (0, i);
}



void trimTrailing (string &s)
{
	size_t i = s. size ();
	while (i)
		if (isSpace (s [i - 1]))
			i--;
		else
			break;
  s. erase (i);
}



void trimLeading (string &s,
                  char c)
{
  size_t i = 0;
  for (; i < s. size () && s [i] == c; i++)
    ;
  s. erase (0, i);
}



void trimTrailing (string &s,
                   char c)
{
	size_t i = s. size ();
	while (i)
		if (s [i - 1] == c)
			i--;
		else
			break;
  s. erase (i);
}



size_t containsWord (const string& hay,
                     const string& needle)
{
  size_t pos = string::npos;
  for (;;)
  {
    pos = hay. find (needle, pos == string::npos ? 0 : pos + 1);
    if (pos == string::npos)
      break;
    const size_t pos1 = pos + needle. size ();
    if (   (pos  == 0            || ! isLetter (hay [pos - 1]))
        && (pos1 == hay. size () || ! isLetter (hay [pos1]))
       )
      break;
  }
  return pos;  
}



void replace (string &s,
              char from,
              char to)
{ 
  for (char& c : s)
	  if (c == from)
	  	c = to;
}



void replace (string &s,
              const string &fromChars,
              char to)
{
  for (char& c : s)
	  if (charInSet (c, fromChars))
	  	c = to;
}



void replaceStr (string &s,
                 const string &from,
                 const string &to)
{
  ASSERT (! from. empty ());
  
  if (from == to)
    return;
  
  const bool inside = contains (to, from);

  size_t start = 0;    
  for (;;)
  {
    const size_t pos = s. find (from, start);
    if (pos == string::npos)
      break;
    s. replace (pos, from. size (), to);
    start = pos + (inside ? to. size () : 0);
  }
}
  
  
  
string strQuote (const string &s,
                 char quote)
{ 
	string s1 (s);
	replaceStr (s1, "\\", "\\\\");
	const char slashQuote [] = {'\\', quote, '\0'};
	replaceStr (s1, string (1, quote), slashQuote);
	return quote + s1 + quote; 
}



string to_c (const string &s)
{
  string r;
  for (const char c : s)
    if (c == '\n')
      r += "\\n";
    else
    {
      if (charInSet (c, "\"\'\\"))
        r += '\\';
      r += c;
    }   
  return r;
}



void collapseSpace (string &s)
{
  string s1;
  do 
  {
    s1 = s;
    replaceStr (s, "  ", " ");
  }
  while (s1 != s);
}
               


  
string str2streamWord (const string &s,
                       size_t wordNum)
{
  istringstream iss (s);
  string word;
  FOR (size_t, i, wordNum + 1)
    if (iss. eof ())
    	return string ();
    else
      iss >> word;
  return word;
}
                


string str2sql (const string &s)
{
	string r = "'";
  for (const char c : s)
	{
	  if (c == '\'')
	  	r += "'";
	  r += c;
	}
	r += "'";
	
	return r;
}



string sql2escaped (const string &s)
{
  string s1;
  for (const char c : s)
  {
    if (charInSet (c, "[]*_%\\"))
      s1 += '\\';
    s1 += c;
  }
  
  return s1;
}



string findSplit (string &s,
                  char c)
{
	const size_t pos = s. find (c);
	if (pos == string::npos)
	{
		const string s1 (s);
		s. clear ();
		return s1;
	}
	const string before (s. substr (0, pos));
	s. erase (0, pos + 1);
	return before;
}



string rfindSplit (string &s,
                   char c)
{
	const size_t pos = s. rfind (c);
	if (pos == string::npos)
	{
		const string s1 (s);
		s. clear ();
		return s1;
	}
	const string after (s. substr (pos + 1));
	s. erase (pos);
	return after;
}



void reverse (string &s)
{
  FFOR (size_t, i, s. size () / 2)
    swap (s [i], s [s. size () - 1 - i]);
}



List<string> str2list (const string &s,
                       char c) 
{
	List<string> res;
	string s1 (s);
	while (! s1. empty ())
	  res << findSplit (s1, c);
	return res;
}



string list2str (const List<string> &strList,
                 const string &sep) 
{
	string s;
	for (const string& e : strList)
	{
		if (! s. empty ())
			s += sep;
	  s += e;
	}
	return s;
}



bool fileExists (const string &fName)
{
  ifstream f (fName. c_str ());
  bool ok = f. good ();
#ifndef _MSC_VER
  if (ok)
    ok = ! directoryExists (fName);
#endif

  return ok;
}



#ifndef _MSC_VER
bool directoryExists (const string &dirName)
{
  DIR* dir = opendir (dirName. c_str ());
  const bool yes = (bool) (dir);
  if (yes)
  {
    if (closedir (dir))
      ERROR;
  }
  return yes;
}
#endif



//

size_t strMonth2num (const string& month)
{
  if (isDigit (month [0]))
  {
    const int m = str2<int> (month);
    ASSERT (m >= 1);
    ASSERT (m <= 12);
    return (size_t) (m - 1);
  }
  
	size_t i = NO_INDEX;
	     if (month == "Jan")  i = 0;
  else if (month == "Feb")  i = 1;
  else if (month == "Mar")  i = 2;
  else if (month == "Apr")  i = 3;
  else if (month == "May")  i = 4;
  else if (month == "Jun")  i = 5;
  else if (month == "Jul")  i = 6;
  else if (month == "Aug")  i = 7;
  else if (month == "Sep")  i = 8;
  else if (month == "Oct")  i = 9;
  else if (month == "Nov")  i = 10;
  else if (month == "Dec")  i = 11;
  else 
  	ERROR;
  	
  return i;
}




// istream

bool getChar (istream &is,
              char &c)
{
  ASSERT (is. good ());

  const int i = is. get ();
  c = EOF;
  if (is. eof ())
  {
    ASSERT (i == c);
    return false;
  }
  ASSERT (i >= 0 && i <= 255);
  c = static_cast<char> (i);

  return true;
}



void skipLine (istream &is)
{
#if 1
  char c;
  while (getChar (is, c) && c != '\n')  // UNIX
    ;
#else
	char c = '\0';
	while (! is. eof () && c != '\n')  // UNIX
  {	
  	ASSERT (is. good ());
	  c = (char) is. get ();
	}
#endif
}



void readLine (istream &is,
               string &s)
{
  s. clear ();
  for (;;)
  {	
    char c;
    if (! getChar (is, c))
      break;
	  if (c == '\n')  // UNIX
	  	break;
	  s += c;
	}
}



string getToken (istream &is,
                 const string &skip,
                 const string &delimeters)
{
  // Skipping skip
  for (;;)
  {
    IMPLY (! is. eof (), is. good ());
    const char c = (char) is. get ();
    if (is. eof ())
      return string ();
    ASSERT (c);
    if (charInSet (c, skip))
      continue;
    if (charInSet (c, delimeters))
      return string (1, c);
    is. unget ();
    break;
  }  
  
  string token;
  for (;;)
  {
    IMPLY (! is. eof (), is. good ());
    const char c = (char) is. get ();
    if (   is. eof () 
        || charInSet (c, skip)
        || charInSet (c, delimeters)
       )
    {
      is. unget ();
      break;
    }
    ASSERT (c);
    token += c;
  }
  ASSERT (! token. empty ());

  return token;
}

 


// Rand

const long Rand::max_ = 2147483647;



ulong Rand::get (ulong max)
{
  ASSERT (max > 0);
  ASSERT (max <= (ulong) max_);
  
  const long up = max_ - (max_ % (long) max);  // Multiple of max
  do run ();
    while (seed >= up);  // P(stay in cycle) < 0.5
  return (ulong) (seed % (long) max);
}



void Rand::qc () const
{
  if (! qc_on)
    return;
  ASSERT (seed > 0);
  ASSERT (seed < max_);
}



void Rand::run ()
{
	qc ();
  
  const long ia = 16807;
  const long iq = 127773;
  const long ir = 2836;

  const long k = seed / iq;
  seed = ia * (seed - k * iq) - ir * k;
  if (seed < 0)
  	seed += max_;
}




// Threads

size_t Threads::threadsToStart = 0;
	
	


// verbose

namespace
{
	int verbose_ = 0;
}


bool verbose (int inc)
{ 
	return Verbose::enabled () ? (verbose_ + inc > 0) : false;
}


Verbose::Verbose (int verbose_arg)
: verbose_old (verbose_)
{
	if (enabled ())
	  verbose_ = verbose_arg;
}


Verbose::Verbose ()
: verbose_old (verbose_)
{
	if (enabled ())
		verbose_++;
}


Verbose::~Verbose ()
{
	if (enabled ())
		verbose_ = verbose_old;
}



Unverbose::Unverbose ()
{
	if (Verbose::enabled ())
	  verbose_--;
}


Unverbose::~Unverbose ()
{
	if (Verbose::enabled ())
	  verbose_++;
}




// Threads

Threads::Threads (size_t threadsToStart_arg)
{ 
  if (! isMainThread ())
	  throw logic_error ("Threads are started not from the main thread");
  if (! empty ())
	  throw logic_error ("Previous threads have not finished");

	threadsToStart = threadsToStart_arg;
	if (threadsToStart >= threads_max)
		throw logic_error ("Too many threads to start");
		
	threads. reserve (threadsToStart);
	
	if (verbose (1))
    cerr << "# Threads started: " << threadsToStart + 1 << endl;
}	



Threads::~Threads ()
{ 
  {
    Progress prog (threads. size () + 1);
    const string step ("consecutive threads finished");
    prog (step);
    for (auto& t : threads) 
    { 
      t. join ();
      prog (step);
  	}
  }
	  
	threads. clear ();
	threadsToStart = 0;
}




// Root

void Root::saveFile (const string &fName) const
{
	if (fName. empty ())
		return;
  
  OFStream f (fName);
  saveText (f);
}




// Named

Named::Named (const string& name_arg)
: name (name_arg) 
{
#ifndef NDEBUG
  if (! goodName (name))
    ERROR_MSG ("Bad name: " + strQuote (name_arg));
#endif
}



void Named::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
    
  ASSERT (goodName (name));
}




// StringVector

StringVector::StringVector (const string &fName,
                            size_t reserve_size)
{
	reserve (reserve_size);
	searchSorted = true;
	ifstream f (fName);
	string s;
	string prev;
	while (f >> s)
	{
	  *this << s;
	  if (s < prev)
	  	searchSorted = false;
	  prev = s;
	}
}




// DisjointCluster

void DisjointCluster::merge (DisjointCluster &other)
// Trying to keep rank's small
{ 
  if (this == & other)
    return;
  DisjointCluster* root1 =        getDisjointCluster ();
  DisjointCluster* root2 = other. getDisjointCluster ();
	if (root1->rankDC > root2->rankDC)
		root2->parentDC = root1;
	else
	{
		root1->parentDC = root2;
		if (root1->rankDC == root2->rankDC)
			root2->rankDC ++;
	}
}



DisjointCluster* DisjointCluster::getDisjointCluster ()
{ 
  ASSERT (parentDC);
	if (this == parentDC)
	  return this;
	parentDC = parentDC->getDisjointCluster ();
	return parentDC;
}




// Progress

size_t Progress::beingUsed = 0;



void Progress::report () const
{ 
  cerr << '\r';
#ifndef _MSC_VER
  cerr << "\033[2K";
#endif
  cerr << n; 
	if (n_max)
		cerr << " / " << n_max;
	if (! step. empty ())
		cerr << ' ' << step;
	cerr << ' ';
	if (! Threads::empty () && ! contains (step, "thread"))
	  cerr << "(main thread) ";
}




// Input

Input::Input (const string &fName,
	            size_t bufSize,
	            uint displayPeriod)
: buf (new char [bufSize])
, ifs (fName)
, is (& ifs)
, prog (0, displayPeriod)  
{ 
  if (! ifs. good ())
    throw runtime_error ("Cannot open: " + strQuote (fName));
  if (! ifs. rdbuf () -> pubsetbuf (buf. get (), (long) bufSize))
  	throw runtime_error ("Cannot allocate buffer to " + strQuote (fName));
}
 


Input::Input (istream &is_arg,
	            uint displayPeriod)
: is (& is_arg)
, prog (0, displayPeriod)  
{ 
  ASSERT (is);
  if (! is->good ())
    throw runtime_error ("Bad input stream");
}
 


void Input::reset ()
{ 
  ASSERT (is == & ifs);

  ifs. seekg (0); 
  ASSERT (ifs. good ());
  prog. reset ();
}




// LineInput

bool LineInput::nextLine ()
{ 
  ASSERT (is);

	if (eof)  
	{ 
		line. clear ();
	  return false;
	}
	
	try 
	{
  	readLine (*is, line);
  	const bool end = line. empty () && is->eof ();

    if (! commentStart. empty ())
    {
      const size_t pos = line. find (commentStart);
      if (pos != string::npos)
        line. erase (pos);
    }
    trimTrailing (line); 

  	eof = is->eof ();
  	lineNum++;
  
  	if (! end)
  		prog ();
  		
  	return ! end;
  }
  catch (const exception &e)
  {
    throw runtime_error ("Reading line " + toString (lineNum) + ":\n" + line + "\n" + e. what ());
  }
}




// ObjectInput

bool ObjectInput::next (Root &row)
{ 
  ASSERT (is);

	row. clear ();

	if (eof)
	  return false;

	row. read (*is);
	lineNum++;

 	eof = is->eof ();
  if (eof)
  {
  	ASSERT (row. empty ());
  	return false;
  }

	prog ();
	
  ASSERT (is->peek () == '\n');
	row. qc ();

  skipLine (*is);

	return true;
}



// CharInput

char CharInput::get ()
{ 
  ASSERT (is);

  ungot = false;
  
	const char c = (char) is->get ();

	eof = is->eof ();
	ASSERT (eof == (c == (char) EOF));

	if (eol)
	{ 
	  lineNum++;
		charNum = 0;
		prog ();
  }
  else
	  charNum++;

	eol = (eof || c == '\n');

	return /*eof ? (char) EOF :*/ c;
}



void CharInput::unget ()
{ 
  ASSERT (is);
	ASSERT (! ungot);
	
	is->unget (); 
	charNum--;
	ungot = true;
}



string CharInput::getLine ()
{
  string s;
  while (! eof)
  {
    const char c = get ();
    if (eol)
      break;
    s += c;
  }
  return s;
}




// Token

void Token::readInput (CharInput &in)
{
	clear ();

  // Skip spaces
	char c = '\0';
	do { c = in. get (); }
	  while (! in. eof && isSpace (c));
	if (in. eof)
	{
	  ASSERT (empty ());
		return;  
  }
		
  charNum = in. charNum;
	if (c == quote)
	{
		type = eText;
		for (;;)
		{ 
			c = in. get (); 
			if (in. eol)
		  //throw CharInput::Error (in, "ending quote");
		    continue;
			if (c == quote)
				break;
			name += c;
		}
	}
	else if (isDigit (c))
	{
		type = eNumber;
		while (! in. eof && isDigit (c))
		{ 
			name += c;
			c = in. get (); 
		}
		if (! in. eof)
			in. unget ();
		num = str2<uint> (name);
		ASSERT (num != numeric_limits<uint>::max ());
	}
	else if (isLetter (c))
	{
		type = eName;
		while (! in. eof && isLetter (c))
		{ 
			name += c;
			c = in. get (); 
		}
		if (! in. eof)
			in. unget ();
	}
	else 
	{
	  type = eDelimiter;
		name = c;
	}	
	ASSERT (! empty ());
}



void Token::qc () const
{
  if (! qc_on)
    return;
  if (! empty ())
  {
  	ASSERT (! name. empty ());
    ASSERT (! contains (name, quote));
  	IMPLY (type != eText, ! contains (name, ' '));
  	IMPLY (type != eNumber, num == noNum);
  	const char c = name [0];
  	switch (type)
  	{ 
  	  case eText:      break;
  		case eNumber:    ASSERT (isDigit (c)); 
  		                 ASSERT (num != noNum);
  		                 break;
  		case eName:      ASSERT (isLetter (c) && ! isDigit (c)); 
  		                 break;
  		case eDelimiter: ASSERT (name. size () == 1); 
  		                 break;
  		default: throw runtime_error ("Unknown type");
  	}
  }
}




// OFStream

void OFStream::open (const string &dirName,
	                   const string &fileName,
	                   const string &extension)
{ 
	ASSERT (! fileName. empty ());
	ASSERT (! is_open ());
	
	string pathName;
	if (! dirName. empty () && ! isDirName (dirName))
	  pathName = dirName + "/";
	pathName += fileName;
	if (! extension. empty ())
		pathName += "." + extension;
		
	ofstream::open (pathName);

	if (! good ())
	  throw runtime_error ("Cannot create file " + strQuote (pathName));
}




// Csv
 
string Csv::getWord ()
{ 
  ASSERT (goodPos ());
  
  size_t start = pos;
  size_t stop = NO_INDEX;
  if (s [pos] == '\"')
  {
    pos++;
    start = pos;
    findChar ('\"');
    ASSERT (s [pos] == '\"');
    stop = pos;
  }
  
  findChar (',');
  if (stop == NO_INDEX)
    stop = pos;
  pos++;
  
  return s. substr (start, stop - start);
}


  
  
void csvLine2vec (const string &line,
                  StringVector &words)
{
  words. clear ();
  Csv csv (line);
  string s;
  while (csv. goodPos ())
  {
    s = csv. getWord ();
    trim (s);
    words << s;
  }
}




// Json

Json::Json (JsonContainer* parent,
            const string& name)
{ 
  ASSERT (parent);
  
  if (const JsonArray* jArray = parent->asJsonArray ())
  {
    ASSERT (name. empty ());  
    var_cast (jArray) -> data << this;
  }
  else if (const JsonMap* jMap = parent->asJsonMap ())
  {
    ASSERT (! name. empty ());  
    ASSERT (! contains (jMap->data, name));
    var_cast (jMap) -> data [name] = this;
  }
}



Token Json::readToken (istream &is)
{
  const string delim ("[]{},:");
  
  string s;
  bool spaces = true;
  bool isC = false;
  while (is)
  {
    char c;
    is. get (c);
    const bool isDelim = charInSet (c, delim);
    if (spaces)
      if (isSpace (c))
        continue;
      else
      {
        if (isDelim)
          return Token (string (1, c), Token::eDelimiter);
        spaces = false;
        isC = c == '\'';
      }
    else
    {
      // Opposite to toStr()
      if (isC)
      {
        if (c == '\'')
          break;
        else
        if (c == '\\')
        {
          is. get (c);
          if (c == 'n')
            c = '\n';
        }
      }
      else
        if (isSpace (c) || isDelim)
        {
          is. unget ();
          break;
        }
    }
      
    s += c;
  }
  
  if (isC)
    s. erase (0, 1);
  else
    ASSERT (! s. empty ());
  
  return Token (s, isC ? Token::eText : Token::eNumber);
}



void Json::parse (istream& is,
                  const Token& firstToken,
                  JsonContainer* parent,
                  const string& name)
{
  string s (firstToken. name);
  strLower (s);
  
  if (firstToken. type == Token::eDelimiter && s == "{")
    new JsonMap (is, parent, name);
  else if (firstToken. type == Token::eDelimiter && s == "[")
    new JsonArray (is, parent, name);
  else if (firstToken. type == Token::eText && s == "null")
    new JsonNull (parent, name);
  else if (firstToken. type == Token::eNumber)
  {
    if (   contains (s, '.')
        || contains (s, 'e')
       )
    {
      uint decimals = 0;
      size_t pos = s. find ('.');
      if (pos != string::npos)
      {
        pos++;
        while (   pos < s. size () 
               && isDigit (s [pos])
              )
        {
          pos++;
          decimals++;
        }
      }
      new JsonDouble (str2<double> (s), decimals, parent, name);
    }
    else
      new JsonInt (str2<int> (s), parent, name); 
  }
  else
    new JsonString (firstToken. name, parent, name);
}



int Json::getInt () const
{ 
  if (! this || asJsonNull ())
    throw runtime_error ("undefined");
  if (const JsonInt* j = asJsonInt ())
    return j->n;
  throw runtime_error ("Not a JsonInt");
}



double Json::getDouble () const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return numeric_limits<double>::quiet_NaN ();
  if (const JsonDouble* j = asJsonDouble ())
    return j->n;
  throw runtime_error ("Not a JsonDouble");
}



string Json::getString () const
{ 
  if (! this || asJsonNull ())
    throw runtime_error ("undefined");
  if (const JsonString* j = asJsonString ())
    return j->s;
  throw runtime_error ("Not a JsonString");
}



const Json* Json::at (const string& name_arg) const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return nullptr;
  if (const JsonMap* j = asJsonMap ())
    return findPtr (j->data, name_arg);
  throw runtime_error ("Not a JsonMap");
}



const Json* Json::at (size_t index) const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return nullptr;
  if (const JsonArray* j = asJsonArray ())
  {
    if (index >= j->data. size ())
      throw runtime_error ("Index out of range");
    else
      return j->data [index];
  }
  throw runtime_error ("Not a JsonArray");
}        



size_t Json::getSize () const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return 0;
  if (const JsonArray* j = asJsonArray ())
    return j->data. size ();
  throw runtime_error ("Not a JsonArray");
}        



// JsonArray


JsonArray::JsonArray (istream& is,
                      JsonContainer* parent,
                      const string& name)
: JsonContainer (parent, name)
{
  bool first = true;
  for (;;)
  {
    Token token (readToken (is));
    if (token. isDelimiter (']'))
      break;
    if (! first)
    {
      ASSERT (token. isDelimiter (','));
      token = readToken (is);
    }
    parse (is, token, this, string());
    first = false;
  }
}



void JsonArray::print (ostream& os) const
{ 
  os << "[";
  bool first = true;
  for (const Json* j : data)
  {
    ASSERT (j);
    if (! first)
      os << ",";
    j->print (os);
    first = false;
  }
  os << "]";
}



// JsonMap

JsonMap::JsonMap ()
{
  ASSERT (! jRoot);
  jRoot = this;
}



JsonMap::JsonMap (const string &fName)
{
  ifstream ifs (fName. c_str ());
  if (! ifs. good ())
    throw runtime_error ("Cannot open: " + strQuote (fName));
  const Token token (readToken (ifs));
  if (! token. isDelimiter ('{'))
    throw runtime_error ("Json file: " + strQuote (fName) + ": should start with '{'");
  parse (ifs);
}



void JsonMap::parse (istream& is)
{
  ASSERT (data. empty ());
  
  bool first = true;
  for (;;)
  {
    Token token (readToken (is));
    if (token. isDelimiter ('}'))
      break;
    if (! first)
    {
      ASSERT (token. isDelimiter (','));
      token = readToken (is);
    }
    ASSERT (! token. name. empty ());
    const Token colon (readToken (is));
    ASSERT (colon. isDelimiter (':'));
    Json::parse (is, readToken (is), this, token. name);
    first = false;
  }
}



JsonMap::~JsonMap ()
{ 
  for (auto& it : data)
    delete it. second;
}



void JsonMap::print (ostream& os) const
{ 
  os << "{";
  bool first = true;
  for (const auto& it : data)
  {
    if (! first)
      os << ",";
    os << toStr (it. first) << ":";
    const Json* j = it. second;
    ASSERT (j);
    j->print (os);
    first = false;
  }
  os << "}";
}




JsonMap* jRoot = nullptr;




// Offset

size_t Offset::size = 0;




//

void exec (const string &cmd,
           const string &logFName)
{
  ASSERT (! cmd. empty ());
  
  if (verbose ())
  	cout << cmd << endl;
  	
	const int status = system (cmd. c_str ());
	if (status)
	{
	  if (! logFName. empty ())
	  {
	    LineInput f (logFName);
	    throw runtime_error (f. getString ());
	  }
		throw runtime_error ("Command failed:\n" + cmd + "\nstatus = " + toString (status));		
	}
}




// FileItemGenerator

FileItemGenerator::FileItemGenerator (size_t progress_displayPeriod,
                                      bool isDir_arg,
                                      const string& fName_arg)
: ItemGenerator (0, progress_displayPeriod)
, isDir (isDir_arg)
, fName( fName_arg)
{ 
  if (isDir)
  { 
    trimSuffix (fName,  "/");
	  char lsfName [4096] = {'\0'};
  #ifdef _MSC_VER
    throw runtime_error ("Windows");  // ??
  #else
    strcpy (lsfName, P_tmpdir);
    strcat (lsfName, "/XXXXXX");
    EXEC_ASSERT (mkstemp (lsfName) != -1);
    ASSERT (lsfName [0]);
    const int res = system (("ls " + fName + " > " + lsfName). c_str ());
  //printf ("res = %d\n", res);
    ASSERT (! res);
    fName = lsfName;
  #endif
  }      
  f. open (fName);
  ASSERT (f. good ()); 
}



bool FileItemGenerator::next (string &item)
{ 
  if (f. eof ())
    return false;
    
	readLine (f, item);
  if (isDir)
  { const size_t pos = item. rfind ('/');
  	if (pos != string::npos)
      item. erase (0, pos + 1);
  }
  trim (item);
  if (item. empty () && f. eof ())
    return false;

  prog (item);
  
	return true;
}    




// Application

Application::Arg::Arg (const string &name_arg,
                       const string &description_arg)
: Named (name_arg)
, description (description_arg)
{
  trim (name);
  ASSERT (! name. empty ());
  ASSERT (! contains (name, ' '));  
  ASSERT (! contains (name, '='));
  ASSERT (! isLeft (name, "-"));
  ASSERT (! description. empty ());
  ASSERT (name != Application::helpS);
  ASSERT (name != Application::versionS);
}



void Application::Arg::qc () const
{
  if (! qc_on)
  	return;
  	
  ASSERT (! name. empty ());
  ASSERT (! description. empty ());
}



void Application::Key::qc () const
{
  if (! qc_on)
  	return;
  Arg::qc ();
  
  IMPLY (app. gnu, name. size () > 1);
  IMPLY (acronym, app. gnu);
  ASSERT (isUpper (var));
  ASSERT (! var. empty () == (app. gnu && ! flag));
  IMPLY (! requiredGroup. empty (), defaultValue. empty ());
}



void Application::Key::saveText (ostream &os) const
{
  os << '-';
	if (app. gnu)
	{
		os << "-" << name;
		if (! flag)
			os << ' ' << var;
	}
	else
	{
	  os << name; 
	  if (! flag)
	  {
	    os << ' ';
	    const bool quoted = defaultValue. empty () || contains (defaultValue, ' ');
	    if (quoted)
	      os << '\"';
	    os << defaultValue;
	    if (quoted)
	      os << '\"';
	  }
	}
}



string Application::Key::getShortHelp () const
{
	string acronymS;
	if (acronym)
	{
		acronymS = string ("-") + acronym;
		if (! flag)
			acronymS += " " + var;
		acronymS += ", ";
	}
  return acronymS + str ();
}



void Application::addKey (const string &name, 
                          const string &argDescription,
                          const string &defaultValue,
                          char acronym,
                          const string &var)
{
  ASSERT (! name. empty ());
  ASSERT (! contains (name2arg, name));
  if (acronym && contains (char2arg, acronym))
  	throw logic_error ("Duplicate option " + strQuote (string (1, acronym)));
  keys << Key (*this, name, argDescription, defaultValue, acronym, var);
  name2arg [name] = & keys. back ();
  if (acronym)
  	char2arg [acronym] = & keys. back ();
}



void Application::addFlag (const string &name,
                           const string &argDescription,
                           char acronym)
{
  ASSERT (! name. empty ());
  ASSERT (! contains (name2arg, name));
  if (acronym && contains (char2arg, acronym))
  	throw logic_error ("Duplicate option " + strQuote (string (1, acronym)));
  keys << Key (*this, name, argDescription, acronym);
  name2arg [name] = & keys. back ();
  if (acronym)
  	char2arg [acronym] = & keys. back ();
}



void Application::addPositional (const string &name,
                                 const string &argDescription)
{
  ASSERT (! contains (name2arg, name));
  IMPLY (gnu, isUpper (name));
  positionals << Positional (name, argDescription);
  name2arg [name] = & positionals. back ();
}



Application::Key* Application::getKey (const string &keyName) const
{
  if (! contains (name2arg, keyName))
    errorExitStr ("Unknown key: " + strQuote (keyName) + "\n\n" + getInstruction ());
    
  Key* key = nullptr;
  if (const Arg* arg = findPtr (name2arg, keyName))  	
    key = var_cast (arg->asKey ());
  if (! key)
    errorExitStr (strQuote (keyName) + " is not a key\n\n" + getInstruction ());
    
  return key;
}



void Application::setPositional (List<Positional>::iterator &posIt,
	                               const string &value)
{
  if (posIt == positionals. end ())
  {
  	if (isLeft (value, "-"))
      errorExitStr (strQuote (value) + " is not a valid option\n\n" + getInstruction ());
  	else
      errorExitStr (strQuote (value) + " cannot be a positional parameter\n\n" + getInstruction ());
  }
  (*posIt). value = value;
  posIt++;
}



void Application::setRequiredGroup (const string &keyName,
                                    const string &requiredGroup)
{ 
  ASSERT (! requiredGroup. empty ());
	getKey (keyName) -> requiredGroup = requiredGroup; 
}



void Application::qc () const
{
  if (! qc_on)
    return;
  
  ASSERT (! description. empty ());
  ASSERT (! version. empty ());
  IMPLY (! needsArg, positionals. empty ());
  
  for (const Positional& p : positionals)
  	p. qc ();
  	
  map<string,size_t> requiredGroups;
  string requiredGroup_prev;
  for (const Key& key : keys)
  {
  	key. qc ();
  	if (! key. requiredGroup. empty () && key. requiredGroup != requiredGroup_prev)
  	{
  		ASSERT (! contains (requiredGroups, key. requiredGroup));
  		requiredGroups [key. requiredGroup] ++;
  	}
  	requiredGroup_prev = key. requiredGroup;
  }
  for (const auto& it : requiredGroups)
  	ASSERT (it. second > 1);
}
	


string Application::getArg (const string &name) const
{
  if (contains (name2arg, name))  
    return name2arg. at (name) -> value;
  throw logic_error ("Parameter " + strQuote (name) + " is not found");
}



bool Application::getFlag (const string &name) const
{
  const string value (getArg (name));
  const Key* key = name2arg. at (name) -> asKey ();
  if (! key || ! key->flag)
    throw logic_error ("Parameter " + strQuote (name) + " is not a flag");
  return value == "true";
}



string Application::key2shortHelp (const string &name) const
{
  if (! contains (name2arg, name))  
    throw logic_error ("Parameter " + strQuote (name) + " is not found");
  const Arg* arg = name2arg. at (name);
  ASSERT (arg);
  if (const Key* key = arg->asKey ())
    return key->getShortHelp ();
  throw logic_error ("Parameter " + strQuote (name) + " is not a key");
}



string Application::getInstruction () const
{
  string instr (description);
  instr += "\n\nUSAGE:   " + programName;

  for (const Positional& p : positionals)
    instr += " " + p. str ();

  string requiredGroup_prev;
  for (const Key& key : keys)
  {
 		if (key. requiredGroup == requiredGroup_prev)
	  	if (key. requiredGroup. empty ())
		    instr += " ";
	  	else
  		  instr += " | ";
  	else
  		if (requiredGroup_prev. empty ())
  			instr += " (";
  		else 
  		{
 			  instr += ") ";
  			if (! key. requiredGroup. empty ())
  				instr += "(";
  		}
	  if (key. requiredGroup. empty ())
	  	instr += "[";
    instr += key. str ();
	  if (key. requiredGroup. empty ())
	  	instr += "]";
  	requiredGroup_prev = key. requiredGroup;
  }
	if (! requiredGroup_prev. empty ())
		instr += ")";
  
//if (! contains (name2arg, "help"))
    instr += "\nHELP:    " + programName + " " + ifS (gnu, "-") + "-" + helpS;
//if (! contains (name2arg, "version"))
	  instr += "\nVERSION: " + programName + " " + ifS (gnu, "-") + "-" + versionS;

  return instr;
}



string Application::getHelp () const
{
  string instr (getInstruction ());

  const string par ("\n    ");

  if (! positionals. empty ())
  {
	  instr += "\n\nOBLIGATORY PARAMETERS:";
	  for (const Positional& p : positionals)
	    instr += "\n" + p. str () + par + p. description;
	}

  if (! keys. empty ())
  {
	  instr += "\n\nOPTIONAL PARAMETERS:";
	  for (const Key& key : keys)
	  {
	    instr += "\n" + key. getShortHelp () + par + key. description;
	    if (gnu && ! key. defaultValue. empty ())
	    	instr += par + "Default: " + key. defaultValue;
	  }
	}
  
  return instr;
}



int Application::run (int argc, 
                      const char* argv []) 
{
  ASSERT (programArgs. empty ());
  

  addDefaultArgs ();
  
  
	try
  { 
  	// programArgs
    for (int i = 0; i < argc; i++)  
    {
    	// gnu: -abc --> -a -b -c ??
    	string key (argv [i]);
    	string value;
    	bool valueP = false;
    	if (   i 
    		  && ! key. empty () 
    		  && key [0] == '-'
    		 )
    	{
	    	const size_t eqPos = key. find ('=');
	    	if (eqPos != string::npos)
	    	{
	    		value = key. substr (eqPos + 1);
	    		key = key. substr (0, eqPos);
	    		valueP = true;
	    	}
	    }
      programArgs. push_back (key);
      if (valueP)
	      programArgs. push_back (value);
    }
    ASSERT (! programArgs. empty ());
	  
      
    // positionals, keys
    {
	    bool first = true;
		  List<Positional>::iterator posIt = positionals. begin ();
	    Key* key = nullptr;
	    Set<const Key*> keysRead;
	    for (string s : programArgs)
	    {
	      if (first)
	      {
	        programName = rfindSplit (s, fileSlash);
	        ASSERT (! programName. empty ());
	      }
	      else if (key)
      	{
          ASSERT (! key->flag);
          key->value = s;
          key = nullptr;
        }
        else if (! s. empty () && s [0] == '-')
        {

          const string s1 (s. substr (1));
          if (s1. empty ())
            errorExitStr ("Dash with no key\n\n" + getInstruction ());

          string name;
          const char c = s1 [0];  // Valid if name.empty()
          if (gnu)
          	if (c == '-')
          	{
          		name = s1. substr (1);
		          if (name. empty ())
		            errorExitStr ("Dashes with no key\n\n" + getInstruction ());
          	}
          	else
          	{
          		if (s1. size () != 1) 
                errorExitStr ("Single character expected: " + strQuote (s1) + "\n\n" + getInstruction ());
            }
          else
          	name = s1;

          if (name == helpS /*&& ! contains (name2arg, helpS)*/)
          {
            cout << getHelp () << endl;
            return 0;
          }
          if (name == versionS /*&& ! contains (name2arg, versionS)*/)
          {
            cout << version << endl;
            return 0;
          }
          
          // key
          if (name. empty ())
          {
          	ASSERT (gnu);
	          key = var_cast (findPtr (char2arg, c));
	        }
          else
          	if (const Arg* arg = findPtr (name2arg, name))					    
					    key = var_cast (arg->asKey ());
					    
					if (key)
          {
	          if (keysRead. contains (key))
	            errorExitStr ("Parameter " + strQuote (key->name) + " is used more than once");
	          else
	            keysRead << key;
	            
	          if (key->flag)
	          {
	            key->value = "true";
	            key = nullptr;
	          }
	        }
	        else
	        	setPositional (posIt, s);
        }
        else
        	setPositional (posIt, s);

	      first = false;
	    }
      if (key)
        errorExitStr ("Key with no value: " + key->name + "\n\n" + getInstruction ());


	    if (programArgs. size () == 1 && (! positionals. empty () || needsArg))
	    {
	      cout << getInstruction () << endl;
	      return 1;
	    }
	    

	    if (posIt != positionals. end ())
	      errorExitStr ("Too few positional parameters\n\n" + getInstruction ());
	  }
    
    
    {
	    map<string,StringVector> requiredGroup2names;
	    map<string,StringVector> requiredGroup2given;
	    for (const Key& key : keys)
	    	if (! key. requiredGroup. empty ())
	    	{
	    		if (! key. value. empty ())
	    		{
	    			const StringVector& given = requiredGroup2given [key. requiredGroup];
		    		if (! given. empty ())
		    			throw runtime_error ("Parameter " + strQuote (key. name) + " is conflicting with " + given. toString (", "));
		    		requiredGroup2given [key. requiredGroup] << strQuote (key. name);
		    	}
		    	requiredGroup2names [key. requiredGroup] << strQuote (key. name);
	    	}
	    for (const auto& it : requiredGroup2names)
	    	if (requiredGroup2given [it. first]. empty ())
	    		throw runtime_error ("One of required parameters is missing: " + it. second. toString (", "));
	  }
  
  
    string logFName;
    string jsonFName;
    if (gnu)
    {
	  	if (getFlag ("debug"))
	  		qc_on = true;
    }
    else
    {
	    logFName = getArg ("log");
	  	ASSERT (! logPtr);
	    if (! logFName. empty ())
	  		logPtr = new ofstream (logFName, ios_base::app);
	  
	  	if (getFlag ("qc"))
	  		qc_on = true;
	
	  	if (getFlag ("noprogress"))
	  		Progress::disable ();
	  	if (getFlag ("profile"))
	  		Chronometer::enabled = true;
	  			
	  	seed_global = str2<ulong> (getArg ("seed"));
	  	if (! seed_global)
	  		throw runtime_error ("Seed cannot be 0");
	
	  	threads_max = str2<size_t> (getArg ("threads"));
	  	if (! threads_max)
	  		throw runtime_error ("Number of threads cannot be 0");
	
	  	if (threads_max > 1 && Chronometer::enabled)
	  		throw runtime_error ("Cannot profile with threads");
	  
	  	jsonFName = getArg ("json");
	  	ASSERT (! jRoot);
	  	if (! jsonFName. empty ())
	  	{
	  		new JsonMap ();
	  	  ASSERT (jRoot);
	  	}
	  	
	  	sigpipe = getFlag ("sigpipe");
	  }
  
  	const Verbose vrb (gnu ? 0 : str2<int> (getArg ("verbose")));
	  	
  
		qc ();
		initEnvironment ();
  	body ();

  
  	if (! jsonFName. empty ())
  	{
  	  ASSERT (jRoot);
  		OFStream f (jsonFName);
      jRoot->print (f);
      delete jRoot;
      jRoot = nullptr;
    }
  
  	if (! logFName. empty ())
  	{
  	  delete logPtr;
  	  logPtr = nullptr;
    }
	}
	catch (const std::exception &e) 
	{ 
	  errorExit (e. what ());
  }


  return 0;
}




// ShellApplication

ShellApplication::~ShellApplication ()
{
	if (! qc_on && ! tmp. empty ())
	  exec ("rm -f " + tmp + "*");  
}



void ShellApplication::initEnvironment () 
{
  ASSERT (tmp. empty ());
  
  // tmp
  if (useTmp)
  {
    char templateS [] = {'/', 't', 'm', 'p', '/', 'X', 'X', 'X', 'X', 'X', 'X', '\0'}; 
  #if 0
    tmp = tmpnam (NULL);
  #else
    mkstemp (templateS);
    tmp = templateS;
  #endif
  	if (tmp. empty ())
  		throw runtime_error ("Cannot generate a temporary file");
  	if (qc_on)
  		cout << tmp << endl;  
  }

  // execDir
	execDir = getProgramDirName ();
	if (execDir. empty ())
		execDir = which (programName);
	ASSERT (isRight (execDir, "/"));
}



string ShellApplication::which (const string &progName) const
{
	if (tmp. empty ())
	  throw runtime_error ("Temporary file is needed");
	
	try { exec ("which " + progName + " 1> " + tmp + ".src 2> /dev/null"); }
	  catch (const runtime_error &)
	    { return string (); }
	LineInput li (tmp + ".src");
	const string s (li. getString ());
	return getDirName (s);
}

	

void ShellApplication::findProg (const string &progName) const
{
	ASSERT (! progName. empty ());
	ASSERT (! contains (progName, '/'));
	ASSERT (isRight (execDir, "/"));
	
	string dir;
	if (! find (prog2dir, progName, dir))
	{
		dir = fileExists (execDir + progName)
		        ? execDir
		        : which (progName);
	  if (dir. empty ())
	  {
	  	cout << "AMRFinder binary " << shellQuote (progName) << " is not found." << endl;
		  cout << "Please make sure that " << shellQuote (progName) << " is in the same directory as " + shellQuote (Common_sp::programName) + " or is in your $PATH." << endl;
	  	exit (1);
	  }	
	  prog2dir [progName] = dir;
	}
	  
	ASSERT (isRight (dir, "/"));
}



string ShellApplication::fullProg (const string &progName) const
{
	string dir;
	EXEC_ASSERT (find (prog2dir, progName, dir));
	ASSERT (isRight (dir, "/"));
	return dir + progName + " ";
}





}


