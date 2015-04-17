//main source file main.cpp to use the class "Rectangle"
#include <iostream>
#include "rectangle.h" //only include the header file of class

int main (){

     Rectangle r(2,3,4,5), s(3,2,4,6); // r and s are object of class Rectangle
     Rectangle * t= & s; //t is pointer to object s
     cout<<"t height: "<<t->GetHeight()<<endl;
     int rarea=r.GetArea();
     int sarea=s.GetArea();
     if (r.GetHeight() * r.GetWidth() > s.GetHeight() * s.GetWidth())
        cout<<"r ";
     else cout<<"s ";
     cout<<"has the greater area "<< sarea <<endl;
}
