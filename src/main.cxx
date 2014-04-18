#include <iostream>

#include "needle.h"
#include "Rendering.h"


#define numNodes 10


BendingNeedleModel needle(18.0, 10, 10.0);



Rendering* r;



/*
  compute with method of [1]:
  1. Chentanez N, Alterovitz R, Ritchie D, Cho L, Hauser KK, Goldberg K, et al. 
  Interactive simulation of surgical needle insertion and steering. 
  ACM Trans Graph;28(3):1-10. Available from: http://portal.acm.org/citation.cfm?id=1531394

  result should be in the form A ap = b with

  A = ( M - dF_dx * dt**2 * beta - dF_dv* dt * gamma)
  b = F + dF_dx * dt * v + dF_dx * dt**2 * ( 0.5 - beta) * a + dF_dv * dt * (1-gamma) * a
   ?= F + dF_dx * ( dt * v + dt**2 * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a
   ?= F + dt * dF_dx * ( v + dt * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a

*/





bool ex = false;

double yoffset = 0.05;
double xrot = 0;
void simulate()
{

  double dt = 0.01;

  static int c = 0;


  for(int i = 0; i < 1 && c < 30000; i++)
  {
    c++;
    //simulateExplicit(dt);
    //simulateImplicit(dt);
    double error = needle.simulateImplicitChentanez(dt);
    std::cout<<c<<" : "<<yoffset<<" : "<<error<<std::endl;
    switch(10) 
    {
      case 0:
        yoffset -= 0.01;
        needle.setBasePosition( Vector(0,yoffset,0));
      break;
      case 10:
        xrot += c < 100 ? 0.001 : -0.001;
        needle.setBaseDirection( Vector(1,xrot,0));
      break;

    }
    if(c%99 == 0 && c < 600) {
      //yoffset -= 0.1;//yoffset;
      //needle.setBasePosition( Vector(0,yoffset,0));
      //x[0][1] += 0.1;
      //x[1][1] += 0.1;
    }
  // for(int i = 0; i < numNodes; i++)
  // {
  // if(debugOut) {
  //   std::cout<<i<<":"<<std::endl;
  //   std::cout<<"x: "<<x[i]<<std::endl;
  //   std::cout<<"f: "<<f[i]<<std::endl;
  //   std::cout<<"v: "<<v[i]<<std::endl;
  //   std::cout<<"a: "<<a[i]<<std::endl;
  // }
  // }

  }

  r->update(needle.getX());
  if(ex) exit(0);
 
}

int main(int argi, char** argv)
{

  if(argi>1) ex = true;

  r = new Rendering();

  needle.addLagrangeModifier(0, Vector(0,1,0));
  needle.addLagrangeModifier(1, Vector(0,1,0));
  int last = needle.getX().size()-1;
  //needle.setSpring( last, Vector( last, 0, 0 ), 0.1 );
  //needle.setSpring( last-1, Vector( last-1, 0, 0 ), 0.1 );
  //needle.setSpring( last-2, Vector( last-2, 0, 0 ), 0.1 );

  r->setup();
  r->setCallback( simulate );
  r->update(needle.getX());
  r->run();  
  

  return 0;
}
