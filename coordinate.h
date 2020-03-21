//coordinate.h
#ifndef COORDINATE_H
#define COORDINATE_H
//using namespace std;
class coordinate{
 public:
  double x; // x-Coordinate
  double y; // y-Coordinate
  // default constr.
  coordinate();
  // single valued constr.
  coordinate(double value);
  // full constructor
  coordinate(double ix, double iy);
  void showcoordinate();
  void set(coordinate inp);
  double dot(coordinate secondone);
};

#endif
