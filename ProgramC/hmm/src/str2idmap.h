#ifndef STR2IDMAP_H
#define STR2IDMAP_H

#include <vector>
#include <string>
#include <map>
using namespace std;

class Str2IdMap {
  map<string, unsigned long> _toId;
  vector<string> _toStr;
public:
  string getStr(unsigned long id) {
    return _toStr[id];
  }
  
  unsigned long getId(string str) {
    map<string, unsigned long>::iterator f = _toId.find(str);
    unsigned long id;
    if (f==_toId.end()) {
      id = _toId.size();
      _toId[str] = id;
      _toStr.push_back(str);
      return id;
    }
    else 
      return f->second;
  }
};

#endif
