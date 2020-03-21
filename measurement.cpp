#include "measurement.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

measurement::measurement(double ia, double ib, int in){
  n=in;
  a=ia;
  b=ib;
  times.resize(in,0);
  values.resize(in,0);
  times.at(0)=0;
  for(int i=1;i<in;i++){
    if (b==1.0)
      {
	times.at(i)=a*i;
      }
    else
      { 
	    times.at(i)=a*pow(b,i);
      }
  }
}


void measurement::print(){
  for (int i=0;i<n-1;i++)
    {
      std::cout<<times.at(i)<<" "<<values.at(i)<<std::endl;
    }
  if(values.at(n-1)!=0)
    {std::cout<<times.at(n-1)<<" "<<values.at(n-1)<<std::endl;}
}

void measurement::print(double scalex, double scaley){
 for (int i=0;i<n-1;i++)
    {
      std::cout<<times.at(i)*scalex<<" "<<values.at(i)*scaley<<std::endl;
    }
 if(values.at(n-1)!=0)
   {std::cout<<times.at(n-1)*scalex<<" "<<values.at(n-1)*scaley<<std::endl;}
}

void measurement::print(double scaley){
 for (int i=0;i<n-1;i++)
    {
      std::cout<<times.at(i)<<" "<<values.at(i)*scaley<<std::endl;
    }
   if(values.at(n-1)!=0)
     {std::cout<<times.at(n-1)<<" "<<values.at(n-1)*scaley<<std::endl;}
}

void measurement::fileprint(double scalex, double scaley, std::string fname){
  std::ofstream outputfile {fname};
  for (int i=0;i<n-1;i++)
    {
      outputfile<<times.at(i)*scalex<<" "<<values.at(i)*scaley<<std::endl;
    }
 if(values.at(n-1)!=0)
   {outputfile<<times.at(n-1)*scalex<<" "<<values.at(n-1)*scaley<<std::endl;}
}

void measurement::fileprint(double scaley, std::string fname){
  std::ofstream outputfile {fname};
  for (int i=0;i<n-1;i++)
    {
      outputfile<<times.at(i)<<" "<<values.at(i)*scaley<<std::endl;
    }
   if(values.at(n-1)!=0)
     {outputfile<<times.at(n-1)<<" "<<values.at(n-1)*scaley<<std::endl;}
}


void measurement::derivativefileprint(double scaley, std::string fname){
  std::ofstream outputfile {fname};
  for (int i=0;i<n-2;i++)
    {
      outputfile<<times.at(i)<<" "<<(values.at(i+1)-values.at(i))/(times.at(i+1)-times.at(i))*scaley<<std::endl;
    }
 if(values.at(n-1)!=0)
   {
     outputfile<<times.at(n-2)<<" "<<(values.at(n-1)-values.at(n-2))/(times.at(n-1)-times.at(n-2))*scaley<<std::endl;
   }
}

void measurement::derivativefileprint(double scalex, double scaley, std::string fname){
  std::ofstream outputfile {fname};
  for (int i=0;i<n-2;i++)
    {
      outputfile<<times.at(i)*scalex<<" "<<(values.at(i+1)-values.at(i))/(times.at(i+1)-times.at(i))*scaley/scalex<<std::endl;
    }
 if(values.at(n-1)!=0)
   {
     outputfile<<times.at(n-2)*scalex<<" "<<(values.at(n-1)-values.at(n-2))/(times.at(n-1)-times.at(n-2))*scaley/scalex<<std::endl;
   }
}

void measurement::logderivativefileprint(double scaley, std::string fname){
  std::ofstream outputfile {fname};
  for (int i=0;i<n-2;i++)
    {
      if (times.at(i)!=0)
	{
	outputfile<<times.at(i)<<" "<<(log10(values.at(i+1))-log10(values.at(i)))/(log10(times.at(i+1))-log10(times.at(i)))<<std::endl;
	}
    }
 if(values.at(n-1)!=0)
   {
     outputfile<<times.at(n-2)<<" "<<(log10(values.at(n-1))-log10(values.at(n-2)))/(log10(times.at(n-1))-log10(times.at(n-2)))<<std::endl;
   }
}

void measurement::logderivativefileprint(double scalex, double scaley, std::string fname){
  std::ofstream outputfile {fname};
  for (int i=0;i<n-2;i++)
    {
      if (times.at(i)!=0)
	{
	outputfile<<times.at(i)*scalex<<" "<<(log10(values.at(i+1))-log10(values.at(i)))/(log10(times.at(i+1))-log10(times.at(i)))<<std::endl;
	}
    }
 if(values.at(n-1)!=0)
   {
     outputfile<<times.at(n-2)*scalex<<" "<<(log10(values.at(n-1))-log10(values.at(n-2)))/(log10(times.at(n-1))-log10(times.at(n-2)))<<std::endl;
   }
}


void measurement::everythingfileprint(double scaley, std::string iname){
  fileprint(scaley,iname);
  derivativefileprint(scaley,iname+".d1");
  logderivativefileprint(scaley,iname+".dlog");
}

void measurement::everythingfileprint(double scalex, double scaley, std::string iname){
  fileprint(scalex,scaley,iname);
  derivativefileprint(scalex,scaley,iname+".d1");
  logderivativefileprint(scalex,scaley,iname+".dlog");
}

void measurement::growsize(){
  
  if (b==1)
    {
      times.push_back(a*n);
      values.push_back(0.0);
    }
  else
    {
      times.push_back(a*pow(b,n));
      values.push_back(0.0);
    }

  n++;

}


