#include <iostream>

#include "generatedCode.h"
#include "mathheader.h"
#include "Rendering.h"


#define numNodes 10
  
std::vector<Vector> x(numNodes);
std::vector<Vector> f(numNodes);
std::vector<Vector> v(numNodes);
std::vector<Vector> a(numNodes);

Rendering* r;

void simulateExplicit()
{
  std::cout<<"updatesimulation"<<std::endl;

  for(int i = 1; i < x.size(); i++)
{
  double m = 0.01;

  if(i != x.size()-1 ) 
  {

  f[i] = Vector(
     -needle_df_dx(x[i-1],x[i],x[i+1]),
     -needle_df_dy(x[i-1],x[i],x[i+1]),
     -needle_df_dz(x[i-1],x[i],x[i+1])
  );
  } else {
    //f[i] = Vector(0,0,0);
    f[i] = (x[i-1] - x[i]) * 0.1;
  }
  

  f[i] += Vector(0,-1,0) * (m * 9.81);
  /*

  v[i] = Vector(
     needle_df_dxx(x[0],x[1],x[2]),
     needle_df_dyy(x[0],x[1],x[2]),
     needle_df_dzz(x[0],x[1],x[2])
  );

  a[i] = Vector(
     needle_df_dxxx(x[0],x[1],x[2]),
     needle_df_dyyy(x[0],x[1],x[2]),
     needle_df_dzzz(x[0],x[1],x[2])
  );
  */


  double dt = 0.02;


  // damping
  v[i] *= 0.9;
 
  // explicit update
  a[i] = f[i] / m * 1.0;

  v[i] = v[i] + a[i] * dt;

  x[i] = x[i] + v[i] * dt;// + a[i] * dt * dt;

}

  std::cout<<"f: "<<f[1]<<std::endl;
  std::cout<<"v: "<<v[1]<<std::endl;
  std::cout<<"a: "<<a[1]<<std::endl;
  r->update(x);
  
}

int main()
{

  r = new Rendering();

  for(int i = 0; i < numNodes; i++)
  {
     x[i] = Vector(i, 0.0f,0);
     f[i] = Vector(0,0,0);
     v[i] = Vector(0,0,0);
     a[i] = Vector(0,0,0);
  }

  r->setup();
  r->setCallback( simulateExplicit );
  r->update(x);
  r->run();  
  

  return 0;
}
