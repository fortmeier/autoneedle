#include <iostream>
#include "generatedCode.h"

int main()
{
  double x = df_dx(-1,0,1, 0,1,0, 0,0,0 );
  std::cout<<"dfdx: "<<x<<std::endl;
  return 0;
}
