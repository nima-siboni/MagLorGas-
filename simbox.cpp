#include "simbox.h"
#include <vector>
#include "obstacle.h"
#include <random>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

static constexpr double epsilonL=0.0000000001;
static constexpr double minfinity=-100000000;
static constexpr double infinity=100000000;

using namespace std;
// Here, by assigning type l t r b to the boundaries in the initilizer, we have assumed that the user always want to have a periodic boundary condition;
// If this is not the case one can implement a switch for the type of the boundary condition where
// reflective
// periodic
// none: which can be implemented with boundaries.resize(0);


static constexpr double pi2 = 6.28318530718;    

SimBox::SimBox(double iL){
  Lx=Ly=iL;
  boundaries.resize(4);
  boundaries[0].set(-1.0*Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,Ly/2.0,'l');
  boundaries[1].set(-1.0*Lx/2.0,Ly/2.0, Lx/2.0,Ly/2.0  ,'t');
  boundaries[2].set(Lx/2.0,Ly/2.0, Lx/2.0,-1.0*Ly/2.0,'r');
  boundaries[3].set(Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,-1.0*Ly/2.0,'b');

}

SimBox::SimBox(double iLx,double iLy){

  Lx=iLx;
  Ly=iLy;

  boundaries.resize(4);
  boundaries[0].set(-1.0*Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,Ly/2.0,'l');
  boundaries[1].set(-1.0*Lx/2.0,Ly/2.0, Lx/2.0,Ly/2.0  ,'t');
  boundaries[2].set(Lx/2.0,Ly/2.0, Lx/2.0,-1.0*Ly/2.0,'r');
  boundaries[3].set(Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,-1.0*Ly/2.0,'b');
}


SimBox::SimBox(double iLx,double iLy,int inObstacles,double iwObstacle,double ilObstacle){
  /* in this routine it is important to use .at() to access the members of obstacle list. It might happen that the last obstacle is set such that it needs an image, then we need either to increase the obstacle_list or exit the program by accessing the out of bound memeber via .at()
   */
  Lx=iLx;
  Ly=iLy;
  nObstacles=inObstacles;
  wObstacle=iwObstacle;
  lObstacle=ilObstacle;
  
  boundaries.resize(4);
  boundaries[0].set(-1.0*Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,Ly/2.0,'l');
  boundaries[1].set(-1.0*Lx/2.0,Ly/2.0, Lx/2.0,Ly/2.0  ,'t');
  boundaries[2].set(Lx/2.0,Ly/2.0, Lx/2.0,-1.0*Ly/2.0,'r');
  boundaries[3].set(Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,-1.0*Ly/2.0,'b');

  obstacle_list.resize(nObstacles);


  
  // std::random_device rdx;  //Will be used to obtain a seed for the random number engine
  // std::random_device rdy;  //Will be used to obtain a seed for the random number engine
  // std::random_device rdtheta;  //Will be used to obtain a seed for the random number engine
  std::mt19937 genx(1); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 geny(100000); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 gentheta(100); //Standard mersenne_twister_engine seeded with rd()
  double hilfe= sqrt(wObstacle*wObstacle/4.0+lObstacle*lObstacle);
  std::uniform_real_distribution<> disx(-1*Lx/2.0, Lx/2.0);
  std::uniform_real_distribution<> disy(-1*Ly/2.0, Ly/2.0);
  std::uniform_real_distribution<> distheta(0,pi2);
  
  double dx=0;
  double dy=0;
  double dtheta=0;
  Obstacle sampleObstacle(wObstacle,lObstacle,0,0);
  int counter=0;
  while (counter<nObstacles){

    dx=disx(genx);
    dy=disy(geny);
    dtheta=distheta(gentheta);

    if (dx>Lx/2.0-hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx-Lx,dy,dtheta);
      counter++;
    }

    if (dx<-1*Lx/2.0+hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx+Lx,dy,dtheta);
      counter++;
    }

    if (dy>Ly/2.0-hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy-Ly,dtheta);
      counter++;
    }

    if (dy<-1*Ly/2.0+hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy+Ly,dtheta);
      counter++;
    }

    if (dx<=Lx/2.0-hilfe && dx>=-1*Lx/2.0+hilfe && dy<=Ly/2.0-hilfe && dy>=-1*Ly/2.0+hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;
    }
    
  }//end of while
  
}

SimBox::SimBox(double iLx,double iLy,int inObstacles,double iwObstacle,double ilObstacle, int randomseed){

  Lx=iLx;
  Ly=iLy;
  nObstacles=inObstacles;
  wObstacle=iwObstacle;
  lObstacle=ilObstacle;
  
  boundaries.resize(4);
  boundaries[0].set(-1.0*Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,Ly/2.0,'l');
  boundaries[1].set(-1.0*Lx/2.0,Ly/2.0, Lx/2.0,Ly/2.0  ,'t');
  boundaries[2].set(Lx/2.0,Ly/2.0, Lx/2.0,-1.0*Ly/2.0,'r');
  boundaries[3].set(Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,-1.0*Ly/2.0,'b');

  obstacle_list.resize(nObstacles);


  /*creating three random number for the seeds using seed_seq*/
  std::seed_seq seq{61396,67454,72163,82624+randomseed,481841+5*randomseed+1,12000-randomseed};
  std::vector<std::uint32_t> seeds(3);
  seq.generate(seeds.begin(), seeds.end());
  std::mt19937 genx(seeds.at(0)); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 geny(seeds.at(1)); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 gentheta(seeds.at(2)); //Standard mersenne_twister_engine seeded with rd()
  double hilfe= sqrt(wObstacle*wObstacle/4.0+lObstacle*lObstacle);
  std::uniform_real_distribution<> disx(-1*Lx/2.0, Lx/2.0);
  std::uniform_real_distribution<> disy(-1*Ly/2.0, Ly/2.0);
  std::uniform_real_distribution<> distheta(0,pi2);
  
  double dx=0;
  double dy=0;
  double dtheta=0;
  Obstacle sampleObstacle(wObstacle,lObstacle,0,0);


  int counter=0;
  while (counter<nObstacles){

    dx=disx(genx);
    dy=disy(geny);
    dtheta=distheta(gentheta);

    if (dx>Lx/2.0-hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx-Lx,dy,dtheta);
      counter++;
    }

    if (dx<-1*Lx/2.0+hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx+Lx,dy,dtheta);
      counter++;
    }

    if (dy>Ly/2.0-hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy-Ly,dtheta);
      counter++;
    }

    if (dy<-1*Ly/2.0+hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;

      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy+Ly,dtheta);
      counter++;
    }
    if (dx<=Lx/2.0-hilfe && dx>=-1*Lx/2.0+hilfe && dy<=Ly/2.0-hilfe && dy>=-1*Ly/2.0+hilfe){
      obstacle_list.at(counter)=sampleObstacle;
      obstacle_list.at(counter).manipulate(dx,dy,dtheta);
      counter++;
    }
    
  }//end of while
  
}


SimBox::SimBox(std::string configfilename){
  ifstream ist {configfilename};
  Obstacle hilfeObstacle;
  std::vector<double> xdata(12);
  std::vector<double> ydata(12);

  double currentxmax,currentxmin,currentymax,currentymin;
  double globalxmin,globalxmax,globalymin,globalymax;
  globalxmax=minfinity;
  globalymax=minfinity;
  globalxmin=infinity;
  globalymin=infinity;
      
  while (ist>>xdata.at(0)>>ydata.at(0)) //if the file is not finisheed
    {
      
      
      for (int i=1;i<12;i++) //read the rest of the obstacle
	{
	  ist>>xdata.at(i);
	  ist>>ydata.at(i);
	}
      
      // for (int i=0; i<12; i++) //read the rest of the obstacle
      // 	{
      // 	  std::cout<<xdata.at(i)<<" "<<ydata.at(i)<<std::endl;
      // 	}
      // std::cout<<xdata.at(0)<<" "<<ydata.at(0)<<std::endl;
      // std::cout<<" "<<std::endl;
      
      
      currentxmax=*max_element(begin(xdata),end(xdata));
      if (currentxmax>globalxmax) {globalxmax=currentxmax;}
      
      currentxmin=*min_element(begin(xdata),end(xdata));
      if (currentxmin<globalxmin) {globalxmin=currentxmin;}

      currentymax=*max_element(begin(ydata),end(ydata));
      if (currentymax>globalymax) {globalymax=currentymax;}

      currentymin=*min_element(begin(ydata),end(ydata));
      if (currentymin<globalymin) {globalymin=currentymin;}


      for (int edgei=0;edgei<12;edgei++)//set the obstacle as the hilfeobstacleg
	{
	  hilfeObstacle.segs[edgei].set(xdata.at(edgei),ydata.at(edgei),xdata.at((edgei+1)%12),ydata.at((edgei+1)%12));
	}
      obstacle_list.push_back(hilfeObstacle);
      //     obstacle_list.at(obstacle_list.size()-1).print();
      
    }


  // setting the private members: Lx and Ly and nObstacles
  Lx = globalxmax-globalxmin+10*epsilonL;
  cout<<"# Lx = "<<Lx<<endl;
  Ly = globalymax-globalymin+10*epsilonL;
  cout<<"# Ly = "<<Ly<<endl;
  nObstacles=obstacle_list.size();
  cout<<"# nObstacles = "<<nObstacles<<endl;

  // translating the obstacles such that they are spread from -Lx/2,Lx/2
  for (int obsti; obsti<obstacle_list.size();obsti++)
    {
      obstacle_list.at(obsti).manipulate(-1.0*(globalxmax+globalxmin)/2.0,-1.0*(globalymax+globalymin)/2.0,0.0);
    }

  // dont forget the boundaries
  
  boundaries.resize(4);
  boundaries[0].set(-1.0*Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,Ly/2.0,'l');
  boundaries[1].set(-1.0*Lx/2.0,Ly/2.0, Lx/2.0,Ly/2.0  ,'t');
  boundaries[2].set(Lx/2.0,Ly/2.0, Lx/2.0,-1.0*Ly/2.0,'r');
  boundaries[3].set(Lx/2.0,-1.0*Ly/2.0,-1.0*Lx/2.0,-1.0*Ly/2.0,'b');

  
}

void SimBox::print(std::string iname){
  std::ofstream outputfile;
  outputfile.open (iname, std::ofstream::out | std::ofstream::trunc);
  outputfile.close();
  for (int i=0; i<nObstacles; i++)
    {
      obstacle_list.at(i).print(iname);
    }
}

double SimBox::getLx(){
  return Lx;
}

double SimBox::getLy(){
  return Ly;
}

int SimBox::getnObstacles(){
  return nObstacles;
}
