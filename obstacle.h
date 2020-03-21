#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "segment.h"
#include <vector>
#include <iostream>
class Obstacle {
 public:
  std::vector<segment> segs;
  Obstacle();
  Obstacle(double w, double l, double centerx, double centery);
  void manipulate(double dx,double dy,double dtheta);//this is something
  void print();
  void print(std::string iname);  
  Obstacle & operator= (const Obstacle & other);
 /* friend: */
 /*  std::ostream& operator<<(std::ostream& os, const Obstacle& obst); */
};// end of the class tracer
#endif
