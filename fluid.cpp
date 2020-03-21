#include <vector>
#include "fluid.h"
#include "tracer.h"
#include "simbox.h"
#include <random>
#include <iostream>
using namespace std;

Fluid::Fluid(int intracers, double inputgyrius, SimBox* simulationbox){

  gyrius=inputgyrius;  
  ntracers=intracers;
  std::mt19937 genx(1000000); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 geny(100000); //Standard mersenne_twister_engine seeded with rd()

  std::uniform_real_distribution<> disx(-1*(simulationbox->getLx())/2.0-gyrius, (simulationbox->getLx())/2.0-gyrius);
  std::uniform_real_distribution<> disy(-1*(simulationbox->getLy())/2.0-gyrius, (simulationbox->getLy())/2.0-gyrius);
  double ranx,rany;
  bool inside;
  if (tracers.size()!=0){
    cout<<"What is wrong! The Fluid object should have been set somewhere else!"<<endl;
  }

  tracers.resize(intracers); //setting the size to intracers.

  for (int tracerid=0; tracerid<tracers.size();tracerid++){
    tracers.at(tracerid).setr(gyrius);
    inside=true;
    while (inside==true){
      ranx=disx(genx);
      rany=disy(geny);
      tracers.at(tracerid).setro(ranx,rany);
      /*check if the tracer is inside an obstacle*/
      int obstid=0;
      bool found=false;
      while ( obstid<simulationbox->getnObstacles() && found==false){
	if (tracers.at(tracerid).isinside_experimental(simulationbox->obstacle_list.at(obstid)) == true){
	  found=true;
	}
	obstid++;
      }
      if (found==false)  inside=false;
    } //end of while
  }//end of loop over tracers
  std::cout<<"#  Succesfully finished initialization of tracers"<<endl;
}


Fluid::Fluid(int intracers, double inputgyrius, SimBox* simulationbox,int randomseed){

  gyrius=inputgyrius;
  ntracers=intracers;
  /*creating two random number for the seeds using seed_seq*/
  std::seed_seq seq{1,20,3200,403,5*randomseed+1,12000,73667,9474+randomseed,19151-randomseed};
  std::vector<std::uint32_t> seeds(3);
  seq.generate(seeds.begin(), seeds.end());
  std::mt19937 genx(seeds.at(0)); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 geny(seeds.at(1)); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 gentheta(seeds.at(2)); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> disx(-1*(simulationbox->getLx())/2.0-gyrius, (simulationbox->getLx())/2.0+gyrius);
  std::uniform_real_distribution<> disy(-1*(simulationbox->getLy())/2.0-gyrius, (simulationbox->getLy())/2.0+gyrius);
  std::uniform_real_distribution<> distheta(0,6.283185);// generate random number in [0,2pi)


  double ranx,rany,rantheta;
  bool inside,insidethesimbox;
  if (tracers.size()!=0){
    cout<<"What is wrong! The Fluid object should have been set somewhere else!"<<endl;
  }

  tracers.resize(intracers); //setting the size to intracers.
  int nrtrials=0;
  for (int tracerid=0; tracerid<tracers.size();tracerid++){
    tracers.at(tracerid).setr(gyrius);
    inside=true; /*inside is true if the tracer is inside an obstacle*/
    insidethesimbox=false;
    while (inside==true || insidethesimbox==false){
      ranx=disx(genx);
      rany=disy(geny);
      rantheta=distheta(gentheta);
      nrtrials++;
      tracers.at(tracerid).setro(ranx,rany);
      tracers.at(tracerid).setangle(rantheta);
      /*check if the tracer is inside the simbox */
      insidethesimbox=false;
      if (tracers.at(tracerid).isinsidesimbox(simulationbox)==true){
	insidethesimbox=true;
      }
      
      /*check if the tracer is inside an obstacle*/
      int obstid=0;
      bool found=false; /* found is true if an obstacle is found such that the tracer is inside it*/
      inside=true; /*inside is true if the tracer is inside an obstacle*/
      while ( obstid<simulationbox->getnObstacles()  && found==false){ 
	if (tracers.at(tracerid).isinside_experimental(simulationbox->obstacle_list.at(obstid)) == true){
	  found=true;
	}
	obstid++;
      }
      /*end of check for obstacles*/
      if (found==false)  inside=false;
    } //end of while
  }//end of loop over tracers
  std::cout<<"# Succesfully finished initialization of tracers"<<std::endl;
  std::cout<<"# with total number of "<<nrtrials<<" trials for initiliazing "<<tracers.size()<<" tracers"<<std::endl;
}









int Fluid::nr(){
  return tracers.size();
}

int Fluid::nrengaged(){
  int inthilfe=0;
  for (int i=0;i<tracers.size();i++)
    {
      if (tracers.at(i).engaged==true)
	{
	  inthilfe++;
	}
    }
  return inthilfe;
}
