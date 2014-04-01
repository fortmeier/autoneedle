#include <iostream>

#include "generatedCode.h"
#include "mathheader.h"
#include "Rendering.h"


#define numNodes 10
  
std::vector<Vector> x(numNodes);

Rendering* r;

void simulate()
{
  std::cout<<"updatesimulation"<<std::endl;
  double dx = needle_jacobian_df_dx(x[0],x[1],x[2]);
  double dy = needle_jacobian_df_dy(x[0],x[1],x[2]);
  double dz = needle_jacobian_df_dz(x[0],x[1],x[2]);

  Vector J(dx,dy,dz);

  x[1] -= J * 0.1;

  std::cout<<"jacobian: "<<J<<std::endl;
  r->update(x);
  
}

int main()
{

  r = new Rendering();

  for(int i = 0; i < numNodes; i++)
  {
     x[i] = Vector(i,pow((double)i,2)/50.0f,0);
  }
  x[1][1] = 3;

  r->setup();
  r->setCallback( simulate );
  r->update(x);
  r->run();  
  

  return 0;
}
