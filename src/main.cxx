#include <iostream>
#include "generatedCode.h"
#include "mathheader.h"

int main()
{
  double x = df_dx(
    Vector(-1,0,0),
    Vector(0,1,0),
    Vector(1,0,0)
  );

  std::cout<<"dfdx: "<<x<<std::endl;
  return 0;
}
