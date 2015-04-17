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

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "tables.h"
#include "logprobs.h"

using std::rand;
const double SMOOTHEDZEROCOUNT=-40;

void TwoDTable::clear()
{
  for (TwoDTable::iterator i = begin(); i!=end(); i++) {
    delete i->second;
  }
  map<unsigned long, OneDTable*>::clear();
  for (map<unsigned long, ULSet*>::iterator j = _possibleContexts.begin(); j!=_possibleContexts.end(); j++) {
    delete j->second;
  }
  _possibleContexts.clear();
  _backoff.clear();
}

TwoDTable::~TwoDTable()
{
  for (TwoDTable::iterator it = begin(); it!=end(); it++) {
    delete it->second;
  }
}

ULSet* TwoDTable::getCntx(unsigned long event)
{
  map<unsigned long, ULSet*>::iterator g = _possibleContexts.find(event);
  if (g==_possibleContexts.end()) 
    return 0;
  else
    return g->second;
}

double TwoDTable::get(unsigned long event, unsigned long context)
{
  TwoDTable::iterator f = find(context);
  if (f==end())
    return _backoff.get(event);
  else
    return f->second->get(event);
}

void TwoDTable::add(unsigned long context, unsigned long event, double count)
{
  TwoDTable::iterator f = find(context);
  OneDTable* entry=0;
  if (f==end()) {
    entry = new OneDTable;
    (*this)[context] = entry;
  }
  else
    entry = f->second;
  entry->add(event, count);
  map<unsigned long, ULSet*>::iterator g = _possibleContexts.find(event);
  ULSet* possCntx=0;
  if (g==_possibleContexts.end()) {
    possCntx = new ULSet;
    _possibleContexts[event] = possCntx;
  }
  else
    possCntx = g->second;
  possCntx->insert(context);
}

void TwoDTable::load(istream& file, Str2IdMap& str2id)
{
  string name1, name2;
  double c;
  
  while (file>>name1>>name2>>c) {
    add(str2id.getId(name1), str2id.getId(name2), log(c));
  }
}

void TwoDTable::save(ostream& file, Str2IdMap& str2id)
{
  for (TwoDTable::iterator i = begin(); i!=end(); i++) {
    OneDTable& vals = *i->second;
    for (OneDTable::iterator j = vals.begin(); j!=vals.end(); j++) {
      file << str2id.getStr(i->first) << ' '
	   << str2id.getStr(j->first) << ' '
	   << exp(j->second) << endl;
    }
  }
}

bool TwoDTable::rand2(unsigned long curr, unsigned long& next)
{
  TwoDTable::iterator it = find(curr);
  if (it==end())
    return false;
  return it->second->rand1(next);
}


OneDTable::OneDTable()
{
  _smoothedZeroCount = SMOOTHEDZEROCOUNT;
}

void OneDTable::load(istream& file, Str2IdMap& str2id)
{
  string name;
  double c;
  while (file>>name>>c) {
    add(str2id.getId(name), log(c));
  }
}

void OneDTable::save(ostream& file, Str2IdMap& str2id)
{
  for (OneDTable::iterator i = begin(); i!=end(); i++) {
      file << str2id.getStr(i->first) << ' '
	   << exp(i->second) << endl;
  }
}

double OneDTable::get(unsigned long event)
{
  OneDTable::iterator f = find(event);
  if (f==end())
    return smoothedZeroCount();
  else
    return f->second;
}

void OneDTable::add(unsigned long event, double count)
{
  OneDTable::iterator f = find(event);
  if (f==end())
    (*this)[event] += count;
  else
    f->second = sumLogProb(f->second, count);
}

bool OneDTable::rand1(unsigned long& next)
{
  double p = double (rand())/RAND_MAX;
  double total = 0;
  for (OneDTable::iterator it = begin(); it!=end(); it++) {
    total += exp(it->second);
    if (total >= p) {
      next = it->first;
      return true;
    }
  }
  return false;
}

