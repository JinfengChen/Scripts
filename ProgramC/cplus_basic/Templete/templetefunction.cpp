#include <iostream>
using namespace std;

template <class T> /*declare the template function. 
                     T is the type of data, 
                     means that all type accepted for this function*/
T abs(T x) /*define the function. T before abs is the type of return data by the function
             T before x is the type of input data of the function */
{
  return x< 0? -x:x;
}

int main()
{
    int n= -5;
    double d = -5.5;
    cout<<abs(n)<<endl;
    cout<<abs(d)<<endl;
}

