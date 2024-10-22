// tsv.hpp

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
*   TSV-table
*
*/


#ifndef TSV_HPP
#define TSV_HPP


#include "common.hpp"
using namespace Common_sp;



namespace Common_sp
{



struct Date : Root
{
  enum Format {fmt_Year, fmt_YMD, fmt_None};  // not complete list ??
  short year {0};
  char month {0};
    // 0 .. 12 - 1
  char day {0};
    // 0 .. days[month] - 1  // leap year ??
  
  
  Date () = default;
  explicit Date (short year_arg,
                 char month_arg = 0,
                 char day_arg = 0)
    : year (year_arg)
    , month (month_arg)
    , day (day_arg)
    {}
  static bool isYear (short n)
    { return n > 1000 && n < 2500; }  // PAR
  static bool isMonth  (short n)
    { return between<short> (n, 0, 12); }  
  static bool isDay  (short n)
    { return between<short> (n, 0, 31); }   // Must depend on month ??
  static Date parse (const string &s,
                     Format fmt);
    // Return: !empty() <=> success
  bool empty () const final
    { return    ! year 
             && ! month
             && ! day; 
    }
  void saveText (ostream &os) const final
    { os << std::setfill('0') << std::setw(4) <<       year      << '-' 
         << std::setfill('0') << std::setw(2) << (int) month + 1 << '-' 
         << std::setfill('0') << std::setw(2) << (int) day   + 1; 
    }
  JsonMap* toJson (JsonContainer* parent, 
                   const string& name = noString) const override
    { auto j = new JsonMap (parent, name);
      new JsonInt (year,  j, "year");
      new JsonInt (month, j, "month");
      new JsonInt (day,   j, "day");
      return j;
    }
    
    
  bool operator== (const Date &other) const
    { return    year  == other. year
             && month == other. month
             && day   == other. day;
    }
  bool less (const Date &other,
             bool equal) const;
  bool operator<= (const Date &other) const
    { return less (other, true); }
  bool operator< (const Date &other) const
    { return less (other, false); }
  Date operator- (const Date &other) const;
    // Requires: other <= *this
  bool year_divisible () const
    { return ! month && ! day; }
  bool quarter_divisible () const
    { return ! (month % 3) && ! day; }
  bool month_divisible () const
    { return ! day; }
};
  
  

struct TextTable : Named
// Tab-delimited (tsv) table with a header
// name: file name or empty()
{
  bool pound {false};
    // '#' in the beginning of header
  bool saveHeader {true};
  struct Header : Named
  { 
    size_t len_max {0};
    // Type
    bool numeric {true};
    // Valid if numeric
    bool scientific {false};
    streamsize decimals {0};
    bool null {false};
      // = can be empty()
    static constexpr size_t choices_max {7};  // PAR
    Set<string> choices;
      // size() <= choices_max + 1
    Header () = default;
    explicit Header (const string &name_arg)
      : Named (name_arg)
      {}
    void qc () const override;
    void saveText (ostream& os) const override
      { os         << name 
           << '\t' << len_max 
           << '\t' << (numeric ? ((scientific ? "float" : "int") + string ("(") + to_string (decimals) + ")") : "char") 
           << '\t' << (null ? "null" : "not null"); 
      }
    void saveSql (ostream& os) const 
      { os << name << ' ';
        if (numeric)
        { if (scientific)
            os << "float";
          else
            os << "numeric(" << len_max << ',' << decimals << ")";
        }
        else
          os << "varchar(" << len_max << ")";
        if (! null)
          os << " not null"; 
        os << endl;
      }
  };
  Vector<Header> header;
    // Header::name's are unique
    // size() = number of columns
  Vector<StringVector> rows;
    // StringVector::size() = header.size()
    // Values are trim()'ed
  typedef  size_t  ColNum;
    // no_index <=> no column
  typedef  size_t  RowNum;
    // no_index <=> no row
  static constexpr char aggr_sep {','};  // PAR
    
    
  struct Error : runtime_error
  {
    Error (const TextTable &tab,
           const string &what)
      : runtime_error (what + "\nIn table file: " + tab. name)
      {}
  };
    

  explicit TextTable (const string &tableFName,
                      const string &columnSynonymsFName = noString);
    // Input: tableFName: format: [{'#' <comment> <EOL>}* '#'] <header> <EOL> {<row> <EOL>>}*
    //                    empty lines are skipped
    //        columnSynonymsFName: <syn_format>
    // Rows where number of columns < header size are added empty values
  static constexpr const char* syn_format {"Column synonyms file with the format: {<main synonym> <eol> {<synonym> <eol>}* {<eol>|<eof>}}*"};
  TextTable () = default;
  TextTable (bool pound_arg,
             const Vector<Header> &header_arg)
    : pound (pound_arg)
    , header (header_arg)
    {}
private:
  void setHeader ();
public:
  void qc () const override;
  void saveText (ostream &os) const override;    
        
  
  static bool getDecimals (string s,
                           bool &hasPoint,
                           streamsize &decimals);
    // Return: true => scientific number
  void printHeader (ostream &os) const;
  ColNum col2num_ (const string &columnName) const;
    // Retuirn: no_index <=> no columnName
  ColNum col2num (const string &columnName) const
    { const ColNum i = col2num_ (columnName);
      if (i == no_index)
        throw Error (*this, "Table has no column " + strQuote (columnName));
      return i;
    }
  Vector<ColNum> columns2nums (const StringVector &columns) const
    { Vector<ColNum> nums;  nums. reserve (columns. size ());
      for (const string &s : columns)
        nums << col2num (s);
      return nums;
    }
  bool hasColumn (const string &columnName) const
    { return col2num_ (columnName) != no_index; }
  void duplicateColumn (const string &columnName_from,
                        const string &columnName_to);
  void substitueColumn (string &columnName_from,
                        const string &columnName_to)
    { duplicateColumn (columnName_from, columnName_to);
      columnName_from = columnName_to;
    }
  ColNum findDate (Date::Format &fmt) const;
    // Date column is not empty and has the same format fmt in all rows
    // Return: no_index <=> not found
    // Output: fmt, valid if return != no_index
  bool isKey (ColNum colNum) const;
private:
  int compare (const StringVector& row1,
               const StringVector& row2,
               ColNum column) const;
public:
  void filterColumns (const StringVector &newColumnNames);
    // Input: newColumnNames: in header::name's
    //          can be repeated
    //          ordered
  void sort (const StringVector &by);
  void group (const StringVector &by,
              const StringVector &sum,
              const StringVector &minV,
              const StringVector &maxV,
              const StringVector &aggr);
    // aggr: slow
    // Invokes: filterColumns(by + sum + aggr)
private:
  void merge (RowNum toRowNum,
              RowNum fromRowNum,
              const Vector<ColNum> &sum,
              const Vector<ColNum> &minV,
              const Vector<ColNum> &maxV,
              const Vector<ColNum> &aggr);
public:
  static StringVector aggr2values (const string &aggr)
    { StringVector v (aggr, aggr_sep, true);
      v. sort ();
      v. uniq ();
      return v;
    }
  void colNumsRow2values (const Vector<ColNum> &colNums,
                          RowNum row_num,
                          StringVector &values) const;
    // Output: values
  RowNum find (const Vector<ColNum> &colNums,
               const StringVector &targetValues,
               RowNum rowNum_start) const;
  StringVector col2values (ColNum col) const;


  struct Key
  {
    const Vector<ColNum> colNums;
    unordered_map<StringVector,RowNum,StringVector::Hasher> data;

    Key (const TextTable &tab,
         const StringVector &columns);
         
    RowNum find (const StringVector &values) const
      { const auto& it = data. find (values);
        if (it != data. end ())
          return it->second;
        return no_index;
      }
  };


  struct Index
  {
    const Vector<ColNum> colNums;
    unordered_map<StringVector,Vector<RowNum>,StringVector::Hasher> data;

    Index (const TextTable &tab,
           const StringVector &columns);
         
    const Vector<RowNum>* find (const StringVector &values) const
      { const auto& it = data. find (values);
        if (it == data. end ())
          return nullptr;
        return & it->second;
      }
  };
};


		
struct TsvOut
{
private:
  ostream* os {nullptr};
  unique_ptr<const ONumber> on;
  size_t lines {0};
  size_t fields_max {0};
  size_t fields {0};
public:
  bool usePound {true};
  
  
  explicit TsvOut (ostream* os_arg,
                   streamsize precision = 6,
	                 bool scientific = false)
	  : os (os_arg)
	  , on (os_arg ? new ONumber (*os_arg, precision, scientific) : nullptr)
	  {}
	  // !os_arg <=> disabled
  explicit TsvOut (ostream &os_arg,
                   streamsize precision = 6,
	                 bool scientific = false)
	  : TsvOut (& os_arg, precision, scientific)
	  {}
 ~TsvOut ()
    { if (fields)
        errorExitStr ("TsvOut: unfinished line with " + to_string (fields) + " fields"); 
    }
    
    
  bool live () const
    { return os; }
  bool empty () const
    { return    ! lines 
             && ! fields_max
             && ! fields;
    }
  template <typename T>
    TsvOut& operator<< (const T &field)
      { if (os)
        { if (lines && fields >= fields_max)
            throw runtime_error ("TsvOut: fields_max = " + to_string (fields_max));
          if (fields)  
            *os << '\t'; 
          else if (! lines && usePound)
            *os << '#';
          *os << field; 
          fields++;
        }
        return *this; 
      }    
  void newLn ()
    { if (! os)
        return;
      *os << endl;
      if (! lines)
        fields_max = fields;
      lines++;
      if (fields != fields_max)
        throw runtime_error ("TsvOut: fields_max = " + to_string (fields_max) + ", but fields = " + to_string (fields));
      fields = 0;
    }
};



}



#endif

