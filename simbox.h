#ifndef SIMBOX_H
#define SIMBOX_H
#include <vector>
#include "segment.h"
#include "obstacle.h"
#include <random>
#include <string>
class SimBox {
 private:
  double Lx;
  double Ly;
  double nObstacles;
  double wObstacle;
  double lObstacle;
 public:
  std::vector<segment> boundaries;
  std::vector<Obstacle> obstacle_list;
  /* SimBox(); */
  SimBox(double iL);
  SimBox(double iLx, double iLy);
  SimBox(double iLx, double iLy, int inObstacles,double iwObstacle,double ilObstacle);
  SimBox(double iLx, double iLy, int inObstacles,double iwObstacle,double ilObstacle, int randomseed);
  SimBox(std::string configfilename);
  double getLx();
  double getLy();
  int getnObstacles();
  void print(std::string iname);


  /* void print(); */
};// end of the class tracer
#endif
