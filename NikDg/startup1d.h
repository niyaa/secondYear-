#ifndef STARTUP1D_H
#define STARTUP1D_H

#include "codes1d.h"
#include "math.h"

class startup1d
{
public:
    startup1d();
    startup1d(int NofOder,  MatrixXd &EToV, VectorXd &VX);


    ~startup1d();

    int N,Nfp,Nfaces,Np,K,TotalFaces,Nv;
    double NODETOL;
    VectorXd r,temp1,temp2;
    MatrixXd nx,Dr,V,invV,EToV,LIFT,x,x1,rx,J,EToE,EToF,Fx,Fscale;



};

#endif // STARTUP1D_H
