// tsv.cpp

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
*   TSV table
*
*/


#undef NDEBUG

#include "tsv.hpp"

#include "common.inc"




namespace Common_sp
{
 


// Date

Date Date::parse (const string &s,
                  Format fmt)
{ 
  istringstream iss (s);
  short year = 0;
  short month = 0;
  short day = 0;
  char c1 = '\0';
  char c2 = '\0';
  string tmp;
  switch (fmt)
  {
    case fmt_Year:
      iss >> year >> tmp;
      if (   tmp. empty ()
          && isYear (year)
         )
        return Date (year);
      break;
    case fmt_YMD:
      iss >> year >> c1 >> month >> c2 >> day >> tmp;
      month--;
      day--;
      if (   tmp. empty ()
          && isYear  (year)
          && isMonth (month)
          && isDay   (day)
          && c1 == c2
         )
        return Date (year, (char) month, (char) day);
      break;
    default: throw runtime_error (FUNC "Unknown date format");
  }
  return Date ();
}



bool Date::less (const Date &other,
                 bool equal) const
{
  LESS_PART (*this, other, year);
  LESS_PART (*this, other, month);
  LESS_PART (*this, other, day);
  return equal;
}



Date Date::operator- (const Date &other) const
{ 
  Date d ( short (year  - other. year)
         , char  (month - other. month)
         , char  (day   - other. day)
         );
  // Normalization 
  // day < 0 ??
  while (d. month < 0)
  {
    d. month = char (d. month + 12);
    d. year --;
  }
  return d;
}




// TextTable

void TextTable::Header::qc () const
{
  if (! qc_on)
    return;

  Named::qc ();
    
  QC_IMPLY (scientific, numeric);
  QC_IMPLY (decimals, numeric);
}



void TextTable::Header::saveSql (ostream& os) const 
{ 
  string name_ (name);
  replace (name_, ' ' ,'_');
  replaceStr (name_, "%", "Percent");
  os << name_ << ' ';
  if (numeric)
  { 
    if (scientific)
      os << "float";
    else
      os << "numeric(" << max<size_t> (1, len_max) << ',' << decimals << ")";
  }
  else
    os << "varchar(" << max<size_t> (1, len_max) << ")";
  if (! null)
    os << " not null"; 
  os << '\n';
}


TextTable::TextTable (const string &tableFName,
                      const string &columnSynonymsFName,
                      bool headerP,
                      uint displayPeriod)
: Named (tableFName)
{  
  {
    LineInput f (tableFName, displayPeriod);
    bool dataExists = true;
    // header
    while (f. nextLine ())
    {
      if (verbose ())
        cerr << f. lineNum << endl;
      trimTrailing (f. line);
      if (f. line. empty ())
        continue;
      const bool thisPound = (f. line. front () == '#');
      if (thisPound)
      {
        pound = true;
        f. line. erase (0, 1);
      }
      if (f. line. empty ())
        continue;
      if (header. empty () || thisPound)
      {
        header. clear ();
        StringVector h (f. line, '\t', true);
        for (string& s : h)
          header << std::move (Header (std::move (s)));
      }
      ASSERT (! header. empty ());
      if (headerP)
      {
        if (! thisPound)
        {
          if (! pound)
            dataExists = f. nextLine ();
          break;
        }
      }
      else
      {
        FFOR (size_t, i, header. size ())
          header [i]. name = to_string (i + 1);
        dataExists = true;
        pound = true;
        break;
      }
    }
    if (header. empty ())
      throw Error (*this, "Cannot read the table header");
    // dataExists <=> f.line is valid
    // rows[]
    while (dataExists)
    {
      trimTrailing (f. line);
      if (! f. line. empty ())
      {
        StringVector row (f. line, '\t', true);
        FFOR_START (size_t, i, row. size (), header. size ())
          row << noString;        
        rows << std::move (row);
        ASSERT (row. empty ());
      }
      dataExists = f. nextLine ();
    }
  }
  
  if (! columnSynonymsFName. empty ())
  {
    LineInput colF (columnSynonymsFName);
    string mainSyn;
    while (colF. nextLine ())
    {
      trim (colF. line);
      const string& syn = colF. line;
      if (syn. empty ())
        mainSyn. clear ();
      else
      {
        if (mainSyn. empty ())
          mainSyn = syn;
        else
          if (mainSyn != syn)
          {
            const ColNum i = col2num_ (syn);
            if (i != no_index)
            {
              if (hasColumn (mainSyn))
                throw runtime_error ("Table " + strQuote (name) + ": Column " + strQuote (mainSyn) + " already exists");
              else
                header [i]. name = mainSyn;
            }
          }
      }
    }
  }
  
  setHeader ();
}



void TextTable::setHeader ()
{
  RowNum row_num = 0;
  for (const StringVector& row : rows)
  {
    row_num++;
    if (row. size () != header. size ())
      throw Error (*this, "Row " + to_string (row_num) + " contains " + to_string (row. size ()) + " columns whereas header has " + to_string (header. size ()) + " columns");
    FFOR (RowNum, i, row. size ())
    {
      string field (row [i]);
      trim (field);
      Header& h = header [i];
      if (strNull (field))
      {
        h. null = true;
        continue;
      }
      maximize (h. len_max, field. size ());
      if (h. choices. size () <= Header::choices_max)
        h. choices << field;
      if (! h. numeric)
        continue;
      {
        char* endptr = nullptr;
        strtod (field. c_str (), & endptr);
        if (endptr != field. c_str () + field. size ())
        {
          h. numeric = false;
          h. scientific = false;
          h. decimals = 0;
        }
      }
      if (h. numeric)
      {
        bool hasPoint = false;
        streamsize decimals = 0;
        if (getScientific (field, hasPoint, decimals))
          h. scientific = true;
        maximize<streamsize> (h. decimals, decimals);
      }
    }
  }

  // Header::len_max for numeric
  for (const StringVector& row : rows)
    FFOR (RowNum, i, row. size ())
    {
      Header& h = header [i];
      string field (row [i]);
      trim (field);
      if (strNull (field))
        continue;
      if (! h. numeric)
        continue;
      bool hasPoint = false;
      streamsize decimals = 0;
      getScientific (field, hasPoint, decimals);
      maximize (h. len_max, field. size () + (size_t) (h. decimals - decimals) + (! hasPoint));
    }
}



void TextTable::qc () const
{
  if (! qc_on)
    return;
  if (! name. empty ())
    Named::qc (); 

  {    
    StringVector v;  v. reserve (header. size ());
    FFOR (size_t, i, header. size ())
    {
      const Header& h = header [i];
      try { h. qc (); }
        catch (const exception &e)
        {
          throw runtime_error ("Header column #" + to_string (i + 1) + ": " + e. what ());
        }
      v << h. name;
    }
    v. sort ();
    const size_t i = v. findDuplicate ();
    if (i != no_index)
      throw Error (*this, "Duplicate column name: " + strQuote (v [i]));
  }
  
  FFOR (RowNum, i, rows. size ())
  {
    if (rows [i]. size () != header. size ())
      throw Error (*this, "Row " + to_string (i + 1) + " contains " + to_string (rows [i]. size ()) + " columns whereas table has " + to_string (header. size ()) + " columns");
    for (const string& field : rows [i])
    {
      if (contains (field, '\t'))
        throw Error (*this, "Field " + strQuote (header [i]. name) + " of row " + to_string (i + 1) + " contains a tab character");
      if (contains (field, '\n'))
        throw Error (*this, "Field " + strQuote (header [i]. name) + " of row " + to_string (i + 1) + " contains an EOL character");
    }
  }
}



void TextTable::saveText (ostream &os) const
{ 
  if (saveHeader)
  { 
    if (pound)
      os << '#';
    bool first = true;
    for (const Header& h : header)
    {
      if (! first)
        os << '\t';
      os << h. name;
      first = false;
    }
    os << endl;
  }
  
  for (const StringVector& row : rows)
  {
    save (os, row, '\t');
    os << endl;
  }
}

    
    
void TextTable::printHeader (ostream &os) const
{
  FFOR (size_t, i, header. size ())
  {
    os << i + 1 << '\t';
    header [i]. saveText (os);
    os << endl;
  }
}



TextTable::ColNum TextTable::col2num_ (const string &columnName) const
{ 
  FFOR (size_t, i, header. size ())
    if (header [i]. name == columnName)
      return i;
  return no_index;
}



void TextTable::duplicateColumn (const string &columnName_from,
                                 const string &columnName_to)
{
  ASSERT (! columnName_to. empty ());
  const ColNum from = col2num (columnName_from);
  if (hasColumn (columnName_to))
    throw runtime_error ("Table already has column " + strQuote (columnName_to));
  header << header [from];
  header. back (). name = columnName_to;
  for (StringVector& row : rows)
    row << row [from];
  qc ();
}



TextTable::ColNum TextTable::findDate (Date::Format &fmt) const
{
  FFOR (ColNum, dateCol, header. size ())
  {
    const Header& h = header [dateCol];
    if (   h. null
        || h. scientific
       )
      continue;
    size_t fmt_ = 0;
    while (fmt_ < Date::fmt_None)
    {
      fmt = Date::Format (fmt_);
      bool isDate = true;
      for (const StringVector& row : rows)
        if (Date::parse (row [dateCol], fmt). empty ())
        {
          isDate = false;
          break;
        }
      if (isDate)
        return dateCol;
      fmt_++;
    }
  }
  return no_index;
}



bool TextTable::isKey (ColNum colNum) const
{
  ASSERT (colNum < header. size ());

  const Header& h = header [colNum];
  if (h. null)
    return false;
  if (h. numeric)
    if (   h. scientific
        || h. decimals > 0
       )
      return false;
    
  unordered_set<string> values;  values. rehash (rows. size ());
  for (const StringVector& row : rows)
  {
    ASSERT (! row [colNum]. empty ());
    if (! values. insert (row [colNum]). second)
      return false;
  }
      
  return true;
}



int TextTable::compare (const StringVector& row1,
                        const StringVector& row2,
                        ColNum column) const
{
  const string& s1 = row1 [column];
  const string& s2 = row2 [column];

  if (header [column]. numeric)
  {
    const double a = strNull (s1) ? 0.0 : stod (s1);
    const double b = strNull (s2) ? 0.0 : stod (s2);
    if (a < b)
      return -1;
    if (a > b)
      return 1;
    return 0;
  }
  
  if (s1 < s2)
    return -1;
  if (s1 > s2)
    return 1;

  return 0;
}



void TextTable::filterColumns (const StringVector &newColumnNames)
{
  const Vector<ColNum> colNums (columns2nums (newColumnNames));

  {  
    Vector<Header> newHeader;  newHeader. reserve (colNums. size ());
    for (const ColNum i : colNums)
      newHeader << header [i];
    header = std::move (newHeader);
  }

  for (StringVector& row : rows)
  {
    StringVector newRow;  newRow. reserve (colNums. size ());
    for (const ColNum i : colNums)
      newRow << row [i];
    row = std::move (newRow);
  }
}



void TextTable::sort (const StringVector &by)
{
  const Vector<ColNum> byIndex (columns2nums (by));
  
  const auto lt = [&byIndex,this] (const StringVector &a, const StringVector &b) 
                    { for (const ColNum i : byIndex) 
                        switch (this->compare (a, b, i))
                        { case -1: return true;
                          case  1: return false;
                        }
                      // Tie resolution 
                      FFOR (ColNum, i, header. size ()) 
                        switch (this->compare (a, b, i))
                        { case -1: return true;
                          case  1: return false;
                        }
                      return false;
                    };
                    
  Common_sp::sort (rows, lt);
}



void TextTable::deredundify (const StringVector &equivCols,
                             const CompareInt& equivBetter)
{
  const Vector<ColNum> byIndex (columns2nums (equivCols));
  
  const auto lt = [&byIndex,this] (const StringVector &a, const StringVector &b) 
                    { for (const ColNum i : byIndex) 
                        switch (this->compare (a, b, i))
                        { case -1: return true;
                          case  1: return false;
                        }
                      return false;
                    };
                    
  Common_sp::sort (rows, lt);
    
  Vector<bool> toDelete (rows. size (), false);
  {
    RowNum i = 0;
    while (i < rows. size ())
    {
      const StringVector& row1 = rows [i];
      FFOR_START (RowNum, j, i + 1, rows. size () + 1)  
      {
        if (j == rows. size ())
        {
          i = rows. size ();
          break;
        }
        const StringVector& row2 = rows [j];
        ASSERT (! lt (row2, row1));
        if (lt (row1, row2))
        {
          i = j;
          break;
        }
        ASSERT (i < j);
        // i and j are in the same equivalence class
        FOR_START (RowNum, k, i, j)
          if (! toDelete [k])
          {
            bool stop = false;
            FOR_START (RowNum, l, k + 1, j + 1)
              if (! toDelete [l])
              {
                // Compare rows[k] and rows[l]
                if (equivBetter (& rows [k], & rows [l]) == 1)
                  toDelete [l] = true;
                else if (equivBetter (& rows [l], & rows [k]) == 1)
                {
                  toDelete [k] = true;
                  stop = true;
                  break;
                }
              }
            if (stop)
              break;
          }
      }
    }
  }
  
  if (exists (toDelete))
  {
    Vector<StringVector> rows_new;  rows_new. reserve (rows. size ());
    FFOR (RowNum, i, rows. size ())
      if (! toDelete [i])
        rows_new << std::move (rows [i]);    
    rows = std::move (rows_new);
  }
}



namespace
{
  struct ColumnPartitionQC
  {
    map<string/*colName*/,string/*operation*/> name2oper;
    
    void add (const StringVector &cols,
              const string &oper)
      {
        for (const string& s : cols)
          if (name2oper [s]. empty ())
            name2oper [s] = oper;
          else
            throw runtime_error ("Column " + strQuote (s) + " is used for operations " + strQuote (name2oper [s]) + " and " + strQuote (oper)); 
      }
  };
}



void TextTable::group (const StringVector &by,
                       const StringVector &sum,
                       const StringVector &minV,
                       const StringVector &maxV,
                       const StringVector &aggr)
{
  const Vector<ColNum> byIndex   (columns2nums (by));
  const Vector<ColNum> sumIndex  (columns2nums (sum));
  const Vector<ColNum> minIndex  (columns2nums (minV));
  const Vector<ColNum> maxIndex  (columns2nums (maxV));
  const Vector<ColNum> aggrIndex (columns2nums (aggr));
  
  // QC
  {
    ColumnPartitionQC cp;
    cp. add (by, "group by");
    cp. add (sum, "sum");
    cp. add (minV, "min");
    cp. add (maxV, "max");
    cp. add (aggr, "aggregation");
  }
  for (const string& s : sum)
    if (! header [col2num (s)]. numeric)
      throw runtime_error ("Summation column " + strQuote (s) + " is not numeric");

  sort (by);
  
  RowNum i = 0;  
  FFOR_START (RowNum, j, 1, rows. size ())
  {
    ASSERT (i < j);
    if (rows [i]. same (rows [j], byIndex))
      merge (i, j, sumIndex, minIndex, maxIndex, aggrIndex);
    else
    {
      i++;
      if (i < j)
        rows [i] = std::move (rows [j]);
    }
  }
  if (! rows. empty ())
    i++;

  ASSERT (rows. size () >= i);
  FFOR (RowNum, k, rows. size () - i)
    rows. pop_back ();
    
  StringVector newColumns;
  newColumns << by << sum << minV << maxV << aggr;
  filterColumns (newColumns);
}



void TextTable::merge (RowNum toRowNum,
                       RowNum fromRowNum,
                       const Vector<ColNum> &sum,
                       const Vector<ColNum> &minV,
                       const Vector<ColNum> &maxV,
                       const Vector<ColNum> &aggr) 
{
  ASSERT (toRowNum < fromRowNum);
  
        StringVector& to   = rows [toRowNum];
  const StringVector& from = rows [fromRowNum];

  for (const ColNum i : sum)
  {
    const Header& h = header [i];
    ASSERT (h. numeric);
    ostringstream oss;
    ONumber on (oss, h. decimals, h. scientific);
    const string& s1 = to   [i];
    const string& s2 = from [i];
    const double d1 = strNull (s1) ? 0.0 : stod (s1);
    const double d2 = strNull (s2) ? 0.0 : stod (s2);
    oss << (d1 + d2);
    to [i] = oss. str ();
  }

  for (const ColNum i : minV)
    if (   to [i]. empty ()
        || (   ! from [i]. empty ()
            && compare (to, from, i) == 1
           )
       )
      to [i] = from [i];

  for (const ColNum i : maxV)
    if (   to [i]. empty ()
        || (   ! from [i]. empty ()
            && compare (to, from, i) == -1
           )
       )
      to [i] = from [i];

  for (const ColNum i : aggr)
  {
    if (from [i]. empty ())
      continue;
    if (contains (from [i], aggr_sep))
      throw runtime_error ("Cannot aggregate column " + header [i]. name + " for row " + to_string (fromRowNum + 1) + " because it contains " + strQuote (string (1, aggr_sep)));
    if (to [i]. empty ())
      to [i] = from [i];
    else
      aggregate (to [i], from [i], aggr_sep);
  }
}



void TextTable::colNumsRow2values (const Vector<ColNum> &colNums,
                                   RowNum row_num,
                                   StringVector &values) const
{
  values. clear ();
  values. reserve (colNums. size ());
  const StringVector& row = rows [row_num];    
  FFOR (ColNum, i, colNums. size ())
    values << row [colNums [i]];
}



TextTable::RowNum TextTable::find (const Vector<ColNum> &colNums,
                                   const StringVector &targetValues,
                                   RowNum row_num_start) const
{
  ASSERT (colNums. size () == targetValues. size ());
  ASSERT (row_num_start != no_index);
  StringVector values;
  FOR_START (RowNum, i, row_num_start, rows. size ())
  {
    colNumsRow2values (colNums, i, values);
    if (values == targetValues)
      return i;
  }
  return no_index;
}



StringVector TextTable::col2values (ColNum col) const
{
  QC_ASSERT (col < header. size ());
  
  Set<string> s;
  for (const StringVector& row : rows)
    if (! row [col]. empty ())
      s << row [col];
            
  StringVector vec;  vec. reserve (s. size ());
  insertAll (vec, s);
  
  return vec;
}




// TextTable::Key


TextTable::Key::Key (const TextTable &tab,
                     const StringVector &columns)
: colNums (tab. columns2nums (columns))
{
  data. rehash (tab. rows. size ());
  StringVector values;  
  FFOR (RowNum, i, tab. rows. size ())
  {
    tab. colNumsRow2values (colNums, i, values);
    for (const string& s : values)
      if (s. empty ())
        throw Error (tab, "Empty value in key, in row " + to_string (i + 1));
    if (data. find (values) != data. end ())
      throw Error (tab, "Duplicate key " + values. toString (",") + " for the key on " + columns. toString (","));
    ASSERT (i != no_index);
    data [values] = i;
  }  
}




// TextTable::Index


TextTable::Index::Index (const TextTable &tab,
                         const StringVector &columns)
: colNums (tab. columns2nums (columns))
{
  data. rehash (tab. rows. size ());
  StringVector values;  
  FFOR (RowNum, i, tab. rows. size ())
  {
    tab. colNumsRow2values (colNums, i, values);
    data [values] << i;
  }  
}




}


