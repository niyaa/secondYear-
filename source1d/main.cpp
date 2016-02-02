#include "source1d.h"
#include "../../usr/include/eigen3/Eigen/Dense"
#include <iostream>

#include "startUp1d.h"
int main()
{

  VectorXd v;
  v.setLinSpaced(5,0,1);

//  std::cout<<JacobiP(v,0,0,4);
//  std::cout<<Vandermonde1D(4,v);
//  std::cout<<JacobiGQ(0,0,4);

int K=4;
int N=3;
int Nv=K+1;





startUp1d ob(N,K);




}
