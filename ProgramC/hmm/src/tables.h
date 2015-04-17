/////////////////////////////////////////////////////////////////////////
//Copyright (C) 2003 Dekang Lin, lindek@cs.ualberta.ca
//
//Permission to use, copy, modify, and distribute this software for any
//purpose is hereby granted without fee, provided that the above
//copyright notice appear in all copies and that both that copyright
//notice and this permission notice appear in supporting documentation.
//No representations about the suitability of this software for any
//purpose is made. It is provided "as is" without express or implied
//warranty.
//
/////////////////////////////////////////////////////////////////////////

#ifndef TABLES_H
#define TABLES_H

#include <set>
#include "str2idmap.h"


/** One dimensional table for storing first degree probabilities or
    frequency counts (in log scale).  */
class OneDTable : public map<unsigned long, double> {
  double _smoothedZeroCount;
public:
  OneDTable();
  double smoothedZeroCount() { return _smoothedZeroCount;}
  double get(unsigned long event);
  void add(unsigned long event, double count);
  void load(istream& file, Str2IdMap& str2id);
  void save(ostream& file, Str2IdMap& str2id);

  /** Randomly generate the next value according to the probabilities
      in this table. Return true if the next value is filled. */
  bool rand1(unsigned long& next);
  ~OneDTable() {}
};

typedef set<unsigned long> ULSet; // set of unsigned long integers

/** Two-dimentional table for storing second degree probabilities or
    frequency counts (in log scale).  */
class TwoDTable : public map<unsigned long, OneDTable*> {
  map<unsigned long, ULSet*> _possibleContexts; 
  OneDTable _backoff;
public:
  /** Retrieve the value associated with event and context. */
  double get(unsigned long event, unsigned long context);

  /** Add the count to the  */
  void add(unsigned long event, unsigned long context, double count);

  /** Return the set of contexts such that (event, context) has an
      associated values. */
  ULSet* getCntx(unsigned long event);

  /** fill the table with data in the file. */
  void load(istream& file, Str2IdMap& str2id);

  /** save the contents in the table in the file. */
  void save(ostream& file, Str2IdMap& str2id);

  /** clear the contents of this table */
  void clear(); 

  /** Randomly generate the next value according to the probabilities
      in this table. Return true if the next value is filled. */
  bool rand2(unsigned long curr, unsigned long& next);

  ~TwoDTable();
};

#endif
