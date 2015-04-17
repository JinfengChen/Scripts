#include <Sequence/Crit.hpp>
#include <Sequence/descriptiveStats.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>

using namespace Sequence;
using namespace std;

//use BOOST random number generation schemes
typedef boost::rand48 RNG;
typedef boost::uniform_01< RNG > _Uni;

int main()
{
  //fill a vector<double> with 1000 random values
  RNG generator(static_cast<int>(std::time(0)));
  _Uni uni(generator);
  vector<double> x;
  for(unsigned i = 0 ; i < 100000 ; ++i)
    {
      x.push_back( uni() );
    }

  //sort the list
  sort(x.begin(),x.end());

  //get upper 95% critical value
  std::pair< std::vector<double>::iterator, double > p = 
    upperCrit()(x.begin(),x.end());
  cout << *(p.first) << '\t' << p.second << endl;

  //get lower 5% critical value
  p = lowerCrit()(x.begin(),x.end());
  cout << *(p.first) << '\t' << p.second << endl;

  //examples of descriptive stats
  std::pair<double,double> desc = meanAndVar(x.begin(),x.end());
  cout << desc.first << '\t' << desc.second << endl;
  cout << mean(x.begin(),x.end()) << '\t'
       << variance(x.begin(),x.end()) << endl;
}
