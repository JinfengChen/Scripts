//source file rectangle.cpp of class "Rectangle"
//in this file the member function in rectangle.h should be implemented
#include "rectangle.h"

Rectangle::Rectangle (int x, int y, int h, int w)
{
   xLow=x;
   yLow=y;
   height=h;
   width =w;
}

int Rectangle::GetHeight(){return height;}
int Rectangle::GetWidth(){return width;}
int Rectangle::GetArea(){return height*width;}
