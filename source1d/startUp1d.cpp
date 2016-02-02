#include "startUp1d.h"
#include "iostream"
startUp1d::startUp1d(int N,int K){
NODETOL=1e-10;

    r = JacobiGL(0,0,N);
    Np = N+1;
    Nfp = 1;
    Nfaces =2;

    V = Vandermonde1D(N,r);
    Dr = Dmatrix1D(N,r,V);
    invV = V.inverse();
    MatrixXd Emat(Np,Nfaces*Nfp);

    Emat.setZero(Np,2);
    Emat(0,0)=1.0;Emat(Np-1,1)=1.0;
    LIFT = V*V.transpose()*Emat;
    EToV.setZero(K,2); VX.setZero(K+1);
    EToV = meshGen(0.0,1.0,K,Nv,VX);

    ArrayXi va,vb;
    va = EToV.col(0);vb=EToV.col(1);
    VectorXd v1(K),v2(K);
{int  temp=0;
    for(int i=0;i<K;i++)
    {temp = va(i);
        v1(i)= VX(temp-1);
        temp = vb(i);
        v2(i)= VX(temp-1);}

    }

VectorXd temp1;
temp1.setOnes(N+1,1);
    x =temp1*v1.transpose()+ 0.5*((r.array()+1).matrix())*(v2-v1).transpose();


    J = Dr*x;
    rx = 1.0/J.array();

Fmask.setZero(2);
    Fmask(0)=0;
    Fmask(1)=Np-1;

    std::cout<<"Fx"<<Fmask<<"\n";

    Fx.setZero(2,K);
    Fx.row(0)=x.row(0);
    Fx.row(1)=x.row(Np-1);
            std::cout<<"Fx"<<Fx<<"\n";

            nx=Normals1D(Nfp,Nfaces,K);
Fscale.setZero(2,K);
Fscale.row(0)=J.row(0);
Fscale.row(1)=J.row(Np-1);

Fscale = 1.0/Fscale.array();


EToEF = Connect1DEToE(EToV);


std::cout<<"Fscale "<<"\n"<<Fscale<<"\n";


std::cout<<"Fscale "<<"\n"<<EToEF<<"\n";






}
