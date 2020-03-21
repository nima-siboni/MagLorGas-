#include <iostream>
#include <string>
#include <cmath>
#include "coordinate.h"
#include "tracer.h"
#include "segment.h"
#include "event.h"
#include "simbox.h"
#include "obstacle.h"
#include "fluid.h"
#include "measurement.h"
#include <vector>

using namespace std;
constexpr double pi = 3.1415926535;    
static constexpr double infinity = 100000000; // 10^-12

int main(int argc, char *argv[])
{

  /****** Input Parameters ******/
  coordinate Origin(0.0,0.0);

  string configinputfile=(argv[1]);
  double lifespan=atof(argv[2]);
  int ntracers=atoi(argv[3]);
  double gyrationradius=atof(argv[4]);
  int randomseed=atoi(argv[5]);
  
  

  /****** Creating the simulation box with the obstacles *****/
  //SimBox simulationbox(Lx,Ly,nrobstacles,w,l,randomseed);
  SimBox simulationbox(configinputfile);
  simulationbox.print("obstacles.dat");
  //  cout<<"# nrinputs = "<<argc<<endl<<"# Lx = "<<Lx<<endl<<"# Ly = "<<Ly<<endl<<"# (l,w) = "<<l<<w<<endl<<"# nrObstacle = "<<nrobstacles<<endl<<"# lifespan = "<<lifespan<<endl<<"# ntracers = "<<ntracers<<endl<<"# gyrationradius ="<<gyrationradius<<endl<<"# randomseed ="<<randomseed<<endl;
  
  /****** Creating the tracers ******/
  Fluid fluid(ntracers,gyrationradius,&simulationbox,randomseed);

  /****** setting up the measurements ******/
  measurement msd(0.01,1.1,2);
    
  double t;
  int timestampid;
  tracer* singletracer;
  tracer tr0;
  
  for (int tracerid=0; tracerid<fluid.nr();tracerid++){ // Loop over tracers
    tr0=fluid.tracers.at(tracerid);
    double life=0;
    timestampid=1;
    event lastevent;
    int progressbarcounter=0;
    singletracer=&(fluid.tracers.at(tracerid));
    while (life < lifespan) { /* Loop for proceeding till the life of the tracers goes beyond the lifespan */
      event nextevent;// It is important to be here insidet the loop over tracers
      nextevent=singletracer->find_nextcollision_single_simbox(&simulationbox);

      if (singletracer->engaged==true)
	{
	  t=nextevent.angle*singletracer->getr();
	  if ( life+t>msd.times.at(timestampid) )
	    {  
	      nextevent.angle=(msd.times.at(timestampid)-life)/singletracer->getr();
	      lastevent.seg=nullptr; // This is very important; if not considered the safemove will mess up 
	      singletracer->move(nextevent.angle);
	      life+=nextevent.angle*singletracer->getr();
	      msd.values.at(timestampid)+=singletracer->distance2(&tr0, &simulationbox);
	      timestampid++;
	      if (timestampid>=msd.n) msd.growsize();
	    }
	  else
	    {
	      if (/* singletracer->engaged == true &&*/ lastevent.seg!=nextevent.seg) // The particle is going to collide with another segment.
		{
		  singletracer->move(nextevent.angle);
		  singletracer->collide(nextevent.seg,&simulationbox);
		  lastevent = nextevent;/* lets hope!*/
		}
	      else // if the particle is colliding with the same segment 
		{
		  if (singletracer->safemove(&nextevent)==true)
		    {
		      singletracer->collide(nextevent.seg,&simulationbox);
		    }
		}

	      life+=nextevent.angle*singletracer->getr();
	    }
	}
      else // If the tracer is not engaged, kill it by increasing its life beyond lifespan
 	{
	  life=lifespan+1;
	}
    } // End of while life 

  }// Loop over tracers

  if (fluid.nrengaged()>0){
    msd.everythingfileprint(1.0/fluid.nrengaged(),"msd");
  }else{
    cout << "# No particle is engaged with any obstacle:\n    # - Either increase the number of tracers, or\n    # - Increase the density of obstacles. "<<endl;
  }
  return 0;
}

  
