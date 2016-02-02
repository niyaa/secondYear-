#include <cstdio>
#include <cstdlib>
#include "all.h"
#include "startup1d.h"
int main(){

int i=6;
MatrixXd EToV(4,2);
EToV<<0,1,1,2,2,3,3,4;
VectorXd VX(5);
VX<<0.0,0.5,1.5,2.5,3.0;
startup1d ob(i,EToV,VX);
std::cout<<ob.nx<<std::endl;
std::cout<<ob.Fscale<<std::endl;
//std::cout<<ob.x<<std::endl;

return 0;


}

