#ifndef MEASUREMENT_H
#define MEASUREMENT_H
#include <vector>
#include <string>

class measurement{
 public:
  int n;
  double a;
  double b;
  
  std::vector<double> times; // 
  std::vector<double> values; // y-Coordinate
  // default constr.
  measurement();
  measurement(double ia,double ib,int in);
  void print();
  void print(double scaley);
  void print(double scalex, double scaley);

  void fileprint(double scaley, std::string iname);
  void fileprint(double scalex, double scaley, std::string iname);
  void derivativefileprint(double scalex, double scaley, std::string iname);
  void derivativefileprint(double scaley, std::string iname);
  void logderivativefileprint(double scaley, std::string iname);
  void logderivativefileprint(double scalex, double scaley, std::string iname);
  void everythingfileprint(double scalex, double scaley, std::string iname);
  void everythingfileprint(double scaley, std::string iname);

  void growsize();
};

#endif
