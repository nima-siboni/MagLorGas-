#ifndef TRACER_H
#define TRACER_H
#include "coordinate.h"
#include "segment.h"
#include "simbox.h"
#include "event.h"
class tracer {
 private:
  coordinate ro; //origin of the gyration circle
  double R;  // radius of the gyration circle
  double theta; //the angle of the particle position
  int nx;
  int ny;
 public:
  bool engaged;
  tracer();
  tracer(coordinate iro,double iR, double itheta);
  tracer(double xo, double yo,double iR, double itheta);
  // return the origin
  coordinate getro();
  /* return the gyration radius*/
  double getr();
  double getangle();
  void setro(coordinate value);
  void setro(double xvalue, double yvalue);

  void setr(double value);
  void setangle(double value);
  void show();
  coordinate getpos();
  coordinate getpos_unfolded(SimBox* simulationbox);
  void showxy();
  double find_nextcollision_single_seg(segment seg);
  event find_nextcollision_single_simbox(SimBox* simulationbox);
  double anglefinder(double thetatarget, double currenttheta);
  void move(double someangle);
  bool safemove(event* nextevent);
  void move_virtual_and_showxy(double someangle, double deltat);
  void collide(segment* seg, SimBox* simbox);
  coordinate getvelocity();
  double velocitydot(double vectorx,double vectory);
  void refindangle(coordinate currentposition);
  void move_virtual_and_showxy(double inputangle);
  /* bool isinside(Obstacle iobstacle); //returns true if the tracer is inside the obstacle */
  bool isinside_experimental(Obstacle iobstacle); //returns true if the tracer is inside the obstacle
  bool isinsidesimbox(SimBox* simbox);
  double distance2(tracer* othertracer,SimBox* simulationbox);
};// end of the class tracer
#endif
