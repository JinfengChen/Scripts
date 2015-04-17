//header file rectangle.h of class "Rectangle"

#include <iostream>
using namespace std;
class Rectangle { // define class 
      public:
            Rectangle (int x, int y, int h, int w); // constructor with parameter
            Rectangle ()
            {xLow=0; yLow=0; height=0; width=0;} // constructor without parameter
            ~Rectangle () {}; // destructor, no need to implemented
            int GetHeight (); // member function to return height
            int GetWidth (); // member function to return width
            int GetArea (); // member function to return area
      private:
            int xLow, yLow, height, width;
};
