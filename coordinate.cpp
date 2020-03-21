#include "coordinate.h"
#include <iostream>

using namespace std;
// default constr.
coordinate::coordinate(){
  x=0.0;
  y=0.0;
}
// single valued constr.
coordinate::coordinate(double value){
    x=value;
    y=value;
}
// full constructor
coordinate::coordinate(double ix, double iy){
  x=ix;
  y=iy;
}

void coordinate::showcoordinate(){
  cout << x <<" "<< y << "\n";
}

void coordinate::set(coordinate input){
  x=input.x;
  y=input.y;
}

double coordinate::dot(coordinate theotherone){
  return x*theotherone.x + y*theotherone.y;
}

ostream& operator<<(ostream& os, const coordinate& cor)  
{  
  cout << cor.x <<" "<< cor.y << "\n";
  return os;  
}  
