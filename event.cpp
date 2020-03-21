#include "event.h"
#include <iostream>
#include "segment.h"

static constexpr double minfinity = -100000000; // 10^-12

using namespace std;

event::event(){
  angle=minfinity;
  seg = nullptr;
}

event::event(double iangle, segment *isegment){
  angle=iangle;
  seg = isegment;
}

