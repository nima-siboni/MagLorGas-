#ifndef EVENT_H
#define EVENT_H
#include "segment.h"

class event{
 public:
  double angle; 
  segment* seg;
  event();
  event(double iangle, segment* isegmentID);
};
#endif
