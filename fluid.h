#ifndef FLUID_H
#define FLUID_H
#include <vector>
#include "simbox.h"
#include "tracer.h"
#include <random>
class Fluid {
 private:
    int ntracers;
    double gyrius;

 public:
    std::vector<tracer> tracers;
    Fluid(int ntracers, double igyrius, SimBox* simulationbox);
    Fluid(int ntracers, double igyrius, SimBox* simulationbox, int randomseed);
    int nr();
    int nrengaged();
};// end of /* TODO:  */he class fluid
#endif
