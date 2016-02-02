#ifndef STARTUP1D_H
#define STARTUP1D_H

#include "all.h"
#include "source1d.h"

class startUp1d{
public:
    startUp1d(){}
    startUp1d(int K,int N);
    ~startUp1d(){}
    startUp1d(int NofOder);
    int N,Nfp,Nfaces,Np,K,Nv;
    double NODETOL;
        VectorXd r,VX;
        MatrixXd x,Dr,V,invV,Lift,Emat,LIFT,rx,J,Fx,Fscale,EToEF;
        MatrixXi EToV,nx;
        ArrayXi Fmask;


};

#endif
