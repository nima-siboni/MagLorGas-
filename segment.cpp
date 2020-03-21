#include "coordinate.h"
#include "segment.h"
#include <iostream>
using namespace std;

segment::segment(){
  rb.x=0.0; // r beginning
  rb.y=0.0; 
  re.x=0.0; // r end
  re.y=0.0;
  segtype='o';
}

segment::segment(coordinate irb, coordinate ire){
  rb.set(irb);
  re.set(ire);
  segtype = 'o';
}

segment::segment(double xb, double yb, double xe, double ye){
  rb.x=xb;
  rb.y=yb;
  re.x=xe;
  re.y=ye;
  segtype = 'o';
}

segment::segment(coordinate irb, coordinate ire, char isegtype){
  rb.set(irb);
  re.set(ire);
  segtype = isegtype;
}

segment::segment(double xb, double yb, double xe, double ye, char isegtype){
  rb.x=xb;
  rb.y=yb;
  re.x=xe;
  re.y=ye;
  segtype = isegtype;
}


void segment::set(coordinate irb, coordinate ire){
  rb.set(irb);                                                                                                                         re.set(ire);
  segtype = 'o';
}

void segment::set(double xb, double yb, double xe, double ye){
  rb.x=xb;
  rb.y=yb;
  re.x=xe;
  re.y=ye;
  segtype = 'o';
}


void segment::set(coordinate irb, coordinate ire, char isegtype){
  rb.set(irb);                                                                                                                         re.set(ire);
  segtype = isegtype;
}

void segment::set(double xb, double yb, double xe, double ye, char isegtype){
  rb.x=xb;
  rb.y=yb;
  re.x=xe;
  re.y=ye;
  segtype = isegtype;
}

ostream& operator<<(ostream& os, const segment& iseg)  
{
  //later: why is this not working?
  // segment::cout<<iseg.rb<<std::endl;
  // segment::cout<<iseg.re<<std::endl;

  cout << iseg.rb.x <<" "<< iseg.rb.y << "\n";
  cout << iseg.re.x <<" "<< iseg.re.y << "\n";

  return os;  
}  
