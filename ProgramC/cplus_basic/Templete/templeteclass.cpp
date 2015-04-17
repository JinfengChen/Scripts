#include <iostream>
#include <cstdlib>
using namespace std;

struct Student
{
   int id;
   float gpa;
};


/* define a simple class*/
class Clock
{
   public:
      Clock(int xx, int yy);
      Clock()
      {x=0;y=0;}
      void ShowTime();
      void SetTime(int xx=0, int yy=0);

   private:
      int x;
      int y;
};

Clock::Clock(int xx, int yy)
{
   x = xx;
   y = yy;
}

inline void Clock::ShowTime()
{
   cout<<x<<":"<<y<<endl;
}

void Clock::SetTime(int xx, int yy)
{
   x = xx;
   y = yy;
}

/* declare the class template*/
template <class T>
class Store
{
   private:
      T item;
      int haveValue;
   public:
      Store(void);
      ~Store(){};
      
      T GetElem(void);
      void PutElem(T x);
};


/* implement the class template*/
template <class T>  // constructor function
Store<T>::Store(void):haveValue(0)
{}

template <class T>  /*function that return the item (type T)*/
T Store<T>::GetElem(void)
{
   if(haveValue == 0)
   {  
     cout<< "No item present!"<<endl;
     exit(1);
   } 
   return item;
}

template <class T> /*function that give the data (type T) to item*/
void Store<T>::PutElem(T x) 
{
   haveValue++;
   item=x;
}


//use the class template
int main()
{
    Student g = {1000,23};
    Store<int> S1, S2;
    Store<Student> S3;
    Store<double> D;
    Store<Clock> N;
    Clock B;  
    Clock C(0,1);


    /* put 3 and -7 to S1 and S2, then output using function in class Store*/
    S1.PutElem(3);
    S2.PutElem(-7);
    cout<<S1.GetElem()<<"  "<<S2.GetElem()<<endl; //output data of S1 and S2

    /* put data struct g to S3, then use GetElem() and id in the data struct to output*/
    S3.PutElem(g);
    cout<<"The student id is"<<S3.GetElem().id<<endl;
    
    /* output empty D by GetElem*/
    //cout<<D.GetElem()<<endl;

    /* test class Clock*/
    B.SetTime(1,2);
    B.ShowTime();
    C.ShowTime();
   
    /* use class Clock as data type in Store*/
    N.PutElem(B); //put a project B as value to Store "N"
    N.GetElem().ShowTime(); //get the item "Clock B" and use ShowTime to display.
                            //note that for whatever input type of PutElem, the output of GetElem will return the same thing
                            //that is to say the function in Store are universal for input data type
}




