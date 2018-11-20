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
ostream* logPtr = nullptr;

bool qc_on = false;
size_t threads_max = 1;
ulong seed_global = 1;

bool Chronometer::enabled = false;



namespace 
{

void segmFaultHandler (int /*sig_num*/)
{
  signal (SIGSEGV, SIG_DFL); 
  errorExit ("Segmentation fault", true); 
}

}



bool initCommon ()
{
  MODULE_INIT

  signal (SIGSEGV, segmFaultHandler);  
        
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
	*os << endl
      << "*** ERROR ***" << endl
      << msg << endl
    #ifndef _MSC_VER
	    << "HOSTNAME: " << (hostname ? hostname : "?") << endl
	    << "PWD: " << (pwd ? pwd : "?") << endl
    #endif
	    << "Progam name: " << programName << endl
	    << "Command line:";
	 for (const string& s : programArgs)
	   *os << " " << s;
	 *os << endl;
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
  if (from == to)
    return;
    
  for (;;)
  {
    const size_t pos = s. find (from);
    if (pos == string::npos)
      break;
  #if 0
    s = s. substr (0, pos) + to + s. substr (pos + from. size ());
  #else
    s. replace (pos, from. size (), to);
  #endif
  }
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
    ERROR_MSG ("Bad name: \"" + name_arg + "\"");
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




// Input

Input::Input (const string &fName,
	            size_t bufSize,
	            uint displayPeriod)
: buf (new char [bufSize])
, ifs (fName)
, prog (0, displayPeriod)  
{ 
  if (! ifs. good ())
    throw runtime_error ("Bad file: '" + fName  + "'");
  if (! ifs. rdbuf () -> pubsetbuf (buf. get (), (long) bufSize))
  	throw runtime_error ("Cannot allocate buffer to '" + fName + "'");
}
 



// LineInput

bool LineInput::nextLine ()
{ 
	if (eof)  
	{ 
		line. clear ();
	  return false;
	}
	
	try 
	{
  	readLine (ifs, line);
    trimTrailing (line); 
  	eof = ifs. eof ();
  	lineNum++;
  
  	const bool end = line. empty () && eof;
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
	row. clear ();

	if (eof)
	  return false;

	row. read (ifs);
	lineNum++;

 	eof = ifs. eof ();
  if (eof)
  {
  	ASSERT (row. empty ());
  	return false;
  }

	prog ();
	
  ASSERT (ifs. peek () == '\n');
	row. qc ();

  skipLine (ifs);

	return true;
}



// CharInput

char CharInput::get ()
{ 
  ungot = false;
  
	const char c = (char) ifs. get ();

	eof = ifs. eof ();
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
	ASSERT (! ungot);
	
	ifs. unget (); 
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
				throw CharInput::Error (in, "ending quote");
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
	  ASSERT (type == eDelimiter);
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
	  throw runtime_error ("Cannot create file \"" + pathName + "\"");
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
    const_cast <JsonArray*> (jArray) -> data << this;
  }
  else if (const JsonMap* jMap = parent->asJsonMap ())
  {
    ASSERT (! name. empty ());  
    ASSERT (! contains (jMap->data, name));
    const_cast <JsonMap*> (jMap) -> data [name] = this;
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
    return NAN;
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
    throw runtime_error ("Bad file: '" + fName + "'");
  const Token token (readToken (ifs));
  if (! token. isDelimiter ('{'))
    throw runtime_error ("Json file: '" + fName + "': should start with '{'");
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
  for (auto it : data)
    delete it. second;
}



void JsonMap::print (ostream& os) const
{ 
  os << "{";
  bool first = true;
  for (const auto it : data)
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

void exec (const string &cmd)
{
  ASSERT (! cmd. empty ());
  EXEC_ASSERT (system (cmd. c_str ()) == 0);
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
  
  ASSERT (! description. empty ());
}



void Application::Key::saveText (ostream &os) const
{
  os << "[-" << name; 
  if (! flag)
  {
    os << " ";
    const bool quoted = value. empty () || contains (value, ' ');
    if (quoted)
      os << "\"";
    os << value;
    if (quoted)
      os << "\"";
  }
  os << "]";
}



void Application::addKey (const string &name, 
                          const string &argDescription,
                          const string &defaultValue)
{
  ASSERT (! contains (args, name));
  keys << Key (name, argDescription, defaultValue);
  args [name] = & keys. back ();
}



void Application::addFlag (const string &name,
                           const string &argDescription)
{
  ASSERT (! contains (args, name));
  keys << Key (name, argDescription);
  args [name] = & keys. back ();
}



void Application::addPositional (const string &name,
                                 const string &argDescription)
{
  ASSERT (! contains (args, name));
  positionals << Positional (name, argDescription);
  args [name] = & positionals. back ();
}



string Application::getArg (const string &name) const
{
  if (contains (args, name))  
    return args. at (name) -> value;
  throw runtime_error ("Parameter \"" + name + "\" is not found");
}



bool Application::getFlag (const string &name) const
{
  const string value (getArg (name));
  const Key* key = args. at (name) -> asKey ();
  if (! key || ! key->flag)
    throw runtime_error ("Parameter \"" + name + "\" is not a flag");
  return value == "true";
}



string Application::getInstruction () const
{
  string instr (description);
  instr += "\nUsage: " + programName;
  for (const Positional& p : positionals)
    instr += " " + p. str ();
  for (const Key& key : keys)
    instr += " " + key. str ();
  
  instr += string ("\n") + "Help:  " + programName + " -help|-h";

  return instr;
}



string Application::getHelp () const
{
  string instr (getInstruction ());
  instr += "\nParameters:";
  const string par ("\n  ");
  for (const Positional& p : positionals)
    instr += par + p. str () + ": " + p. description;
  for (const Key& key : keys)
    instr += par + key. str () + ": " + key. description;;
  
  return instr;
}



int Application::run (int argc, 
                      const char* argv []) 
{
  ASSERT (programArgs. empty ());
	try
  { 
    for (int i = 0; i < argc; i++)  
      programArgs. push_back (argv [i]);
    ASSERT (! programArgs. empty ());
  
      
    // positionals, keys
    bool first = true;
    posIt = positionals. begin ();
    Key* key = nullptr;
    Set<string> keysRead;
    for (string s : programArgs)
    {
      if (first)
      {
        programName = rfindSplit (s, fileSlash);
        ASSERT (! programName. empty ());
      }
      else
      {
        if (! s. empty () && s [0] == '-')
        {
          if (key)
            errorExitStr ("Key with no value: " + key->name + "\n" + getInstruction ());
          const string name (s. substr (1));
          if (   name == "help"
          	  || name == "h"
          	 )
          {
            cout << getHelp () << endl;
            return 0;
          }
          if (! contains (args, name))
            errorExitStr ("Unknown key: " + name + "\n" + getInstruction ());
          key = const_cast <Key*> (args [name] -> asKey ());
          if (! key)
            errorExitStr (name + " is not a key\n" + getInstruction ());
          if (keysRead. contains (name))
            errorExitStr ("Parameter \"" + name + "\" is used more than once");
          else
            keysRead << name;
          if (key->flag)
          {
            key->value = "true";
            key = nullptr;
          }
        }
        else
          if (key)
          {
            ASSERT (! key->flag);
            key->value = s;
            key = nullptr;
          }
          else
          {
            if (posIt == positionals. end ())
              errorExitStr ("Too many positional arguments\n" + getInstruction ());
            (*posIt). value = s;
            posIt++;
          }
      }
      first = false;
    }
  
  
    const string logFName = getArg ("log");
  	ASSERT (! logPtr);
    if (! logFName. empty ())
  		logPtr = new ofstream (logFName, ios_base::app);
  
  	if (getFlag ("qc"))
  		qc_on = true;

  	Verbose vrb (str2<int> (getArg ("verbose")));
  	
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
  
  	const string jsonFName = getArg ("json");
  	ASSERT (! jRoot);
  	if (! jsonFName. empty ())
  	{
  		new JsonMap ();
  	  ASSERT (jRoot);
  	}
  
  
    if (programArgs. size () == 1 && (! positionals. empty () || needsArg))
    {
      cout << getInstruction () << endl;
      return 1;
    }
    
    if (posIt != positionals. end ())
      errorExitStr ("Too few positional arguments\n" + getInstruction ());


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




}


