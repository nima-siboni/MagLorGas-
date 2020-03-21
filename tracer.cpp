#include "coordinate.h"
#include "segment.h"
#include "tracer.h"
#include "event.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

static constexpr double pi = 3.1415926535;    // define them as static so that every creation of the tracers do not create them
static constexpr double pi2 = 6.28318530718;    
static constexpr double pihalf = 1.5707963268;    
static constexpr double epsilonL = 0.0000001; // 10^-7
static constexpr double infinity = 100000000; // 10^-12
static constexpr double minfinity = -100000000; // 10^-12

tracer::tracer(){
  ro.x = 0.0; //ro : r origin
  ro.y = 0.0;
  R=0.0; // gyration radius
  theta=0.0; // the angle
  engaged = true; //default value should be true
  nx=ny=0;
}

tracer::tracer(coordinate iro,double iR, double itheta){
  ro.x = iro.x;
  ro.y = iro.y;
  R=iR;
  theta=itheta;
  engaged = true; //default value should be true
  nx=ny=0;
}

tracer::tracer(double xo, double yo,double iR, double itheta){
  ro.x = xo;
  ro.y = yo;
  R=iR;
  theta=itheta;
  engaged = true; //default value should be true
  nx=ny=0;
}
// return the origin
coordinate tracer::getro(){
  return ro;
}
/* return the gyration radius*/
double tracer::getr(){ 
  return R;
}

double tracer::getangle(){
  return theta;
}
void tracer::setro(coordinate value){
  ro.x = value.x;
  ro.y = value.y;
}
void tracer::setro(double xvalue, double yvalue){
  ro.x = xvalue;
  ro.y = yvalue;
}

void tracer::setr(double value){
  R=value;
}

void tracer::setangle(double value){
  
  theta = value - floor(value/pi2)*pi2;
  
}

void tracer::show(){
  cout << ro.x<<" "<<ro.y<<" "<<R<<" "<<theta<< "\n";
}

coordinate tracer::getpos(){
  coordinate temp(ro.x+R*cos(theta),ro.y+R*sin(theta));
  return temp;
}

coordinate tracer::getpos_unfolded(SimBox* simulationbox){
  coordinate temp(ro.x+nx*simulationbox->getLx()+R*cos(theta),ro.y+ny*simulationbox->getLy()+R*sin(theta));
  return temp;
}

void tracer::showxy(){
  coordinate temp;
  temp=getpos();
  cout << temp.x <<" "<<temp.y<< "\n";
}

coordinate tracer::getvelocity(){
 coordinate temp(cos(theta+pihalf),sin(theta+pihalf));
 return temp;
}
  
double  tracer::velocitydot(double vectorx, double vectory){
 return vectorx*cos(theta+pihalf)+vectory*sin(theta+pihalf);
}
  
  
double tracer::anglefinder(double thetatarget, double currenttheta)// finds how far should the tracer move to reach the target from current;
{
  return (thetatarget-currenttheta>0) ? thetatarget-currenttheta : pi2 + thetatarget-currenttheta ;
}

double  tracer::find_nextcollision_single_seg(segment seg){
  double a = pow( seg.rb.x - seg.re.x , 2 ) + pow( seg.rb.y - seg.re.y , 2 );
  double minusb = -2 * ( (seg.re.x - seg.rb.x)*(seg.rb.x - ro.x ) + (seg.re.y - seg.rb.y)*(seg.rb.y - ro.y ) );
  double c = pow(seg.rb.x-ro.x,2) + pow(seg.rb.y-ro.y,2) - R*R;
  double delta = (pow(minusb,2)-4*a*c);
  coordinate rc1;
  coordinate rc2;
  double theta1,theta2;
  theta1=0.0; //initialize to zero to check later if this collision has really happened.
  theta2=0.0; //initialize to zero to check later if this collision has really happened.
  bool collision1 = false;
  bool collision2 = false;
  double collisionangle1,collisionangle2,collisionangle;
  if ( delta < 0 ) {// return a negative value for no collision; I am only going to accept s in [0,1]
    return  minfinity; 
  }
  else
    {
      double sqdelta = pow(delta,0.5);
      double s1 = (minusb + sqdelta)/(2*a);
      if (s1>0 && s1<1) {
	rc1.x=seg.rb.x+s1*(seg.re.x - seg.rb.x);
        rc1.y=seg.rb.y+s1*(seg.re.y - seg.rb.y);
	double cosv = (rc1.x-ro.x)/R;
	double sinv = (rc1.y-ro.y)/R;
	if (sinv>0){
	  theta1 = acos(cosv);
	}
	else{
	  theta1 = pi2 - acos(cosv);
	}
	collision1 = true;
      } // end of if (s1>0 ...
      
      double s2 = (minusb - sqdelta)/(2*a);
      if (s2>0 && s2<1) { 
	rc2.x=seg.rb.x+s2*(seg.re.x - seg.rb.x);
        rc2.y=seg.rb.y+s2*(seg.re.y - seg.rb.y);
	double cosv = (rc2.x-ro.x)/R;
	double sinv = (rc2.y-ro.y)/R;
	if (sinv>0)
	  {theta2 = acos(cosv);}
	else
	  {theta2 = pi2 - acos(cosv);
	  }
	collision2 = true;
      } // end of if (s2>0 ...
      
      // comment: up to here we have found the two angles of collision; the next step is to find out which of the are happening first.
      if (collision1 == true ) {collisionangle1=anglefinder(theta1,theta);} else {collisionangle1=minfinity;}
      if (collision2 == true ) {collisionangle2=anglefinder(theta2,theta);} else {collisionangle2=minfinity;}

      if ( collision1 == true && collision2 == true){ // if both collisions happen
	if (collisionangle1 > collisionangle2){
	  collisionangle=collisionangle2;
	}else {
	  collisionangle = collisionangle1;
	}
      }// end of if both collision happens
      else //one or zero collision are happening 
	if (collision1==true) collisionangle= collisionangle1; //if the collision one is happening
	else
	  if (collision2==true) collisionangle= collisionangle2; //if the collision two is happening
	  else 
	    collisionangle = minfinity; //if non of them are happening
	
      
      // collisionangle -= epsilon/R;
      return collisionangle;

    }//end of else for delta<0
}

void tracer::move(double inputangle){
  if (inputangle < 0) {cerr<<"something went wrong?!"; } else { 
    theta = theta + inputangle;
    theta = theta - floor(theta/pi2)*pi2;
  }
}


bool tracer::safemove(event *nextevent){

  if (nextevent->angle < 0) {cerr<<"something went wrong?!"; } else {
    bool result=true;    
    double dx=(nextevent->seg->re.x-nextevent->seg->rb.x);
    double dy=(nextevent->seg->re.y-nextevent->seg->rb.y);
    double nx=dy;
    double ny=-1.0*dx;
    double hilfe_before=velocitydot(nx,ny);
    theta = theta + nextevent->angle;
    theta = theta - floor(theta/pi2)*pi2;
    double hilfe_after=velocitydot(nx,ny);
    if (hilfe_before*hilfe_after < 0)
      {result=true;}
    else
      {result=false;}
    return result;
  }
}


void tracer::move_virtual_and_showxy(double inputangle){
  double thetaprime;
  if (inputangle < 0) {cerr<<"fuck off!"; } else { 
    thetaprime = theta + inputangle;
    thetaprime = thetaprime - floor(thetaprime/pi2)*pi2;
  }
  cout<< ro.x+R*cos(thetaprime)<<" "<<ro.y+R*sin(thetaprime)<<endl;

}



void tracer::move_virtual_and_showxy(double inputangle, double deltat){
  double thetaprime=theta;
  int nvirtualsteps=0;
  if (inputangle < 0) 
    {
      cerr<<"fuck off!"; 
    } 
  else 
    { 
      nvirtualsteps= (int)floor(inputangle/deltat);
      for (int i=0; i < nvirtualsteps; i++)
	{
	  thetaprime += deltat;
	  thetaprime = thetaprime - floor(thetaprime/pi2)*pi2;
	  cout<< ro.x+R*cos(thetaprime)<<" "<<ro.y+R*sin(thetaprime)<<endl;
	}
    }
  cout<< ro.x+R*cos(thetaprime)<<" "<<ro.y+R*sin(thetaprime)<<endl;
      
}


void tracer::refindangle(coordinate currentposition){
  double deltax=currentposition.x-ro.x;
  double deltay=currentposition.y-ro.y;
  double cosv = (deltax)/sqrt(deltax*deltax+deltay*deltay);
  double sinv = (deltay)/sqrt(deltax*deltax+deltay*deltay);
  if (sinv>0)
    {theta = acos(cosv);}
  else
    {theta = pi2 - acos(cosv);
    }
  
}

void tracer::collide(segment* seg, SimBox* simbox){ //assigns new center of rotation and the corresponding new angle
  double lx,ly,neworiginx,neworiginy,hilfeQ;
  coordinate position_before_collision;
  switch (seg->segtype){ 

  case 'o':{
    lx = seg->re.x - seg->rb.x;
    ly = seg->re.y - seg->rb.y;
    coordinate q(R*cos(theta),R*sin(theta));
    position_before_collision = getpos();
    hilfeQ = (q.x*lx+q.y*ly)/(lx*lx+ly*ly);
    
    neworiginx = ro.x + 2*hilfeQ*lx;
    neworiginy = ro.y + 2*hilfeQ*ly;
    ro.x = neworiginx;
    ro.y = neworiginy;
    refindangle(position_before_collision);
    break;}

  case 'r':{
    position_before_collision = getpos();
    position_before_collision.x += epsilonL ; // just to put it out the box; because the collison angle includes an 1*epsilon which makes sure that the particles dont go through the obstacles. But here we need it to go beyond the wall. The 1*epsilon is in form of theta; here we need something in form of length that is why it is multiplied by R
    position_before_collision.x -= simbox->getLx() ;
    ro.x -= simbox->getLx()-epsilonL;//moving the 
    //    refindangle(position_before_collision);
    nx++;
    break;}

  case 'l':{
    position_before_collision = getpos();
    position_before_collision.x += simbox->getLx()-epsilonL;
    ro.x += simbox->getLx()-epsilonL;//moving the origin
    //refindangle(position_before_collision);
    nx--;
    break;}

  case 'b':{
    position_before_collision = getpos();
    position_before_collision.y += simbox->getLy()-epsilonL  ; // similar to case 'r'

    ro.y += simbox->getLy()-epsilonL;//moving the 
    //refindangle(position_before_collision);
    ny--;
    break;}

  case 't':{
    position_before_collision = getpos();
    position_before_collision.y -= simbox->getLy()- epsilonL ; // similar to case 'r'
    ro.y -= simbox->getLy()-epsilonL;//moving the 
    //refindangle(position_before_collision);
    ny++;
    break;}

  default:{
    cout<<"What are you doing!\n The type of the segment is not in o,t,b,l,r"<<endl;
    return;}
  }  // end of switch

}



event tracer::find_nextcollision_single_simbox(SimBox* simulationbox){

  double anghilfe=infinity;
  Obstacle* cobst; //current obstacle
  event nextevent;
  double ang;
  if (engaged==true){
    anghilfe=infinity;
    for (int obstacleID=0; obstacleID<simulationbox->obstacle_list.size(); obstacleID++){ /*loop obstacles in the simulation box*/
      cobst=&(simulationbox->obstacle_list.at(obstacleID));//ptr to the current obstacle
      for (int segID=0; segID<cobst->segs.size();segID++){ // going through the segments of the current obstacle
	ang=find_nextcollision_single_seg(cobst->segs.at(segID));
	if (ang>0 && ang<anghilfe) {
	  anghilfe=ang;
	  nextevent.angle=ang;
	  nextevent.seg = &(cobst->segs.at(segID)); //returns the pointer to the segment
	}
      }// end of for segs of the current obstacle
    } // end of the loops for obstacles
  
    for (int segid=0; segid < simulationbox->boundaries.size() ; segid++){
      ang=find_nextcollision_single_seg(simulationbox->boundaries[segid]);
      if (ang>0 && ang<anghilfe) {
	anghilfe=ang;
	nextevent.angle=ang;
	nextevent.seg=&(simulationbox->boundaries[segid]);
      }
    }// end of the loop for segments of the simulation box
  } else {
    nextevent.angle=pi2;
  }

  if ( nextevent.seg == nullptr || nextevent.angle==minfinity ){ // if no collision is found
    engaged = false;
    nextevent.angle=pi2;
  }
  return nextevent; //engaged or not it delivers the nextevent.
}



// bool  tracer::isinside(Obstacle iobstacle){
//   // first each segment is presented as f(x,y)=0 where f(x,y)=y-a*x-c
//   // then f(xtracer,ytracer) is evaluated for segments in front of each other and the sign is checked
//   double xt = getpos().x; // position of the tracer
//   double yt = getpos().y; // position of the tracer

  
//   double dx,dy,f0,f6,f1,f11,f3,f9,f2,f4,a,c; 

  
//   dx=iobstacle.segs.at(0).re.x-iobstacle.segs.at(0).rb.x;
//   dy=iobstacle.segs.at(0).re.y-iobstacle.segs.at(0).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(0).rb.y*(dx)-iobstacle.segs.at(0).rb.x*(dy))/(dx);
//       f0=yt-a*xt-c;
//     }
//   else
//     {
//       f0=xt-iobstacle.segs.at(0).rb.x;
//     }


//   dx=iobstacle.segs.at(6).re.x-iobstacle.segs.at(6).rb.x;
//   dy=iobstacle.segs.at(6).re.y-iobstacle.segs.at(6).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(6).rb.y*(dx)-iobstacle.segs.at(6).rb.x*(dy))/(dx);
//       f6=yt-a*xt-c;
//     }
//   else
//     {
//       f6=xt-iobstacle.segs.at(6).rb.x;
//     }


//   dx=iobstacle.segs.at(1).re.x-iobstacle.segs.at(1).rb.x;
//   dy=iobstacle.segs.at(1).re.y-iobstacle.segs.at(1).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(1).rb.y*(dx)-iobstacle.segs.at(1).rb.x*(dy))/(dx);
//       f1=yt-a*xt-c;
//     }
//   else
//     {
//       f1=xt-iobstacle.segs.at(1).rb.x;
//     }





//   dx=iobstacle.segs.at(11).re.x-iobstacle.segs.at(11).rb.x;
//   dy=iobstacle.segs.at(11).re.y-iobstacle.segs.at(11).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(11).rb.y*(dx)-iobstacle.segs.at(11).rb.x*(dy))/(dx);
//       f11=yt-a*xt-c;
//     }
//   else
//     {
//       f11=xt-iobstacle.segs.at(11).rb.x;
//     }


//   dx=iobstacle.segs.at(3).re.x-iobstacle.segs.at(3).rb.x;
//   dy=iobstacle.segs.at(3).re.y-iobstacle.segs.at(3).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(3).rb.y*(dx)-iobstacle.segs.at(3).rb.x*(dy))/(dx);
//       f3=yt-a*xt-c;
//     }
//   else
//     {
//       f3=xt-iobstacle.segs.at(3).rb.x;
//     }



//   dx=iobstacle.segs.at(9).re.x-iobstacle.segs.at(9).rb.x;
//   dy=iobstacle.segs.at(9).re.y-iobstacle.segs.at(9).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(9).rb.y*(dx)-iobstacle.segs.at(9).rb.x*(dy))/(dx);
//       f9=yt-a*xt-c;
//     }
//   else
//     {
//       f9=xt-iobstacle.segs.at(9).rb.x;
//     }


//   dx=iobstacle.segs.at(2).re.x-iobstacle.segs.at(2).rb.x;
//   dy=iobstacle.segs.at(2).re.y-iobstacle.segs.at(2).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(2).rb.y*(dx)-iobstacle.segs.at(2).rb.x*(dy))/(dx);
//       f2=yt-a*xt-c;
//     }
//   else
//     {
//       f2=xt-iobstacle.segs.at(2).rb.x;
//     }


//   dx=iobstacle.segs.at(4).re.x-iobstacle.segs.at(4).rb.x;
//   dy=iobstacle.segs.at(4).re.y-iobstacle.segs.at(4).rb.y;
//   if (dx!=0)
//     {
//       a=(dy)/(dx);
//       c=(iobstacle.segs.at(4).rb.y*(dx)-iobstacle.segs.at(4).rb.x*(dy))/(dx);
//       f4=yt-a*xt-c;
//     }
//   else
//     {
//       f4=xt-iobstacle.segs.at(4).rb.x;
//     }

  
//   if ( ((f0*f6<0)&&(f1*f11<0)) || ((f3*f9<0)&&(f2*f4<0)) )
//     {
//       return true;
//     }
//   else
//     {
//       return false;
//     }
    
// }



bool  tracer::isinside_experimental(Obstacle iobstacle){
  // first each segment is presented as f(x,y)=0 where f(x,y)=y-a*x-c
  // then f(xtracer,ytracer) is evaluated for segments in front of each other and the sign is checked
  double xt = getpos().x; // position of the tracer
  double yt = getpos().y; // position of the tracer

  
  double dx,dy,f1,f7,f0,f2,f10,f4,f11,f9,a,c;




  dx=iobstacle.segs.at(1).re.x-iobstacle.segs.at(1).rb.x;
  dy=iobstacle.segs.at(1).re.y-iobstacle.segs.at(1).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(1).rb.y*(dx)-iobstacle.segs.at(1).rb.x*(dy))/(dx);
      f1=yt-a*xt-c;
    }
  else
    {
      f1=xt-iobstacle.segs.at(1).rb.x;
    }



  
  dx=iobstacle.segs.at(7).re.x-iobstacle.segs.at(7).rb.x;
  dy=iobstacle.segs.at(7).re.y-iobstacle.segs.at(7).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(7).rb.y*(dx)-iobstacle.segs.at(7).rb.x*(dy))/(dx);
      f7=yt-a*xt-c;
    }
  else
    {
      f7=xt-iobstacle.segs.at(7).rb.x;
    }




  
  dx=iobstacle.segs.at(0).re.x-iobstacle.segs.at(0).rb.x;
  dy=iobstacle.segs.at(0).re.y-iobstacle.segs.at(0).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(0).rb.y*(dx)-iobstacle.segs.at(0).rb.x*(dy))/(dx);
      f0=yt-a*xt-c;
    }
  else
    {
      f0=xt-iobstacle.segs.at(0).rb.x;
    }



  
  dx=iobstacle.segs.at(2).re.x-iobstacle.segs.at(2).rb.x;
  dy=iobstacle.segs.at(2).re.y-iobstacle.segs.at(2).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(2).rb.y*(dx)-iobstacle.segs.at(2).rb.x*(dy))/(dx);
      f2=yt-a*xt-c;
    }
  else
    {
      f2=xt-iobstacle.segs.at(2).rb.x;
    }




  
  dx=iobstacle.segs.at(10).re.x-iobstacle.segs.at(10).rb.x;
  dy=iobstacle.segs.at(10).re.y-iobstacle.segs.at(10).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(10).rb.y*(dx)-iobstacle.segs.at(10).rb.x*(dy))/(dx);
      f10=yt-a*xt-c;
    }
  else
    {
      f10=xt-iobstacle.segs.at(10).rb.x;
    }




  
  dx=iobstacle.segs.at(4).re.x-iobstacle.segs.at(4).rb.x;
  dy=iobstacle.segs.at(4).re.y-iobstacle.segs.at(4).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(4).rb.y*(dx)-iobstacle.segs.at(4).rb.x*(dy))/(dx);
      f4=yt-a*xt-c;
    }
  else
    {
      f4=xt-iobstacle.segs.at(4).rb.x;
    }



  
  dx=iobstacle.segs.at(11).re.x-iobstacle.segs.at(11).rb.x;
  dy=iobstacle.segs.at(11).re.y-iobstacle.segs.at(11).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(11).rb.y*(dx)-iobstacle.segs.at(11).rb.x*(dy))/(dx);
      f11=yt-a*xt-c;
    }
  else
    {
      f11=xt-iobstacle.segs.at(11).rb.x;
    }




  
  dx=iobstacle.segs.at(9).re.x-iobstacle.segs.at(9).rb.x;
  dy=iobstacle.segs.at(9).re.y-iobstacle.segs.at(9).rb.y;
  if (dx!=0)
    {
      a=(dy)/(dx);
      c=(iobstacle.segs.at(9).rb.y*(dx)-iobstacle.segs.at(9).rb.x*(dy))/(dx);
      f9=yt-a*xt-c;
    }
  else
    {
      f9=xt-iobstacle.segs.at(9).rb.x;
    }



  
  if ( ((f1*f7<0)&&(f0*f2<0)) || ((f10*f4<0)&&(f11*f9<0)) )
    {
      return true;
    }
  else
    {
      return false;
    }
    
}


bool  tracer::isinsidesimbox(SimBox* simbox){

  double xt = getpos().x; // position of the tracer
  double yt = getpos().y; // position of the tracer

  double Lxmax,Lxmin,Lymax,Lymin;
  Lxmax= 1*(simbox->getLx())/2.0;
  Lxmin=-1*(simbox->getLx())/2.0;
  Lymax= 1*(simbox->getLy())/2.0;
  Lymin=-1*(simbox->getLy())/2.0;
  if ( (xt<Lxmax) && (xt>Lxmin)  &&    (yt<Lymax) && (yt>Lymin)  )
    {
      return true;
    }
  else
    {
      return false;
    }
    
}





double tracer::distance2(tracer* othertracer, SimBox* simulationbox){
  double dx=getpos_unfolded(simulationbox).x-othertracer->getpos_unfolded(simulationbox).x;
  double dy=getpos_unfolded(simulationbox).y-othertracer->getpos_unfolded(simulationbox).y;
  return dx*dx+dy*dy;
}


