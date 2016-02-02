#include "startup1d.h"

startup1d::startup1d()
{

}

startup1d::~startup1d()
{

}

startup1d::startup1d(int NofOder, MatrixXd &EToV,VectorXd &VX){
    N=NofOder;
    r=JacobiGL(0,0,N);
    NODETOL = 1e-10;
    Np=N+1;
    Nfp=1;K=EToV.rows();Nv=K+1;
    Nfaces=2;TotalFaces=Nfaces*K;
    V=Vandermonde1D(N,r);
    Dr=Dmatrix1D(N,r,V);
    MatrixXd Emat(Np,Nfaces*Nfp);
    Emat(0,0)=1.0;Emat(Np-1,1)=1.0;
       int k=0;
    //Codes for finding the boudary points of r
    ArrayXi fmask1(Nfp),fmask2(Nfp+2);
    double l;
   k=0;
   for (int i=0;i<Np;i++)
    {l=fabs(r(i)+1);
       if(l<NODETOL)
          {
       fmask1(k)=i;
       k=k+1;}
    }
   std::cout<<"fmask1"<<fmask1<<"\n";
k=0;
  for (int i=0;i<Np;i++)
    {l=fabs(r(i)-1);
       if(l<NODETOL)
          {
fmask2(k)=i;
 k=k+1;}
   }
  int p=0,g;
  std::cout<<"fmask2"<<fmask2<<"\n";
    k=fmask1.size()+fmask2.size();
    ArrayXi Fmask(k);
         Fmask<<fmask1,fmask2;
for (int i=0;i<k;i++)
{
    g=Fmask(i);
    std::cout<<x.row(g)<<std::endl;
    Fx.row(p)=x.row(g);
Fscale.row(p)=J.row(g);
p=p+1;
}
Fscale = 1.0/Fscale.array();

std::cout<<"Fscale \n"<<Fscale<<std::endl;

//Normas matrix
MatrixXd n;
n.setOnes(1,K);
nx.setZero(Nfp*Nfaces,K);
nx.row(0)=n;nx.row(1)=n;

std::cout<<nx<<std::endl;

   // LIFT MARRIX
LIFT=V*(V.transpose()*Emat);
////nx.setZero(Nfp*Nfaces,K);
VectorXd va,vb;
MatrixXd c;
va=EToV.col(0);
vb=EToV.col(1);

x.setZero(N+1,K);
c.setOnes(N+1,1);

for (int i=0;i<K;i++)
{
    va(i)=VX(va(i));
    vb(i)=VX(vb(i));
}

// JACOBIAN AND RX
J.setZero(N+1,K);
x= c*va.transpose()+0.5*(r+c)*(vb-va).transpose();


J=Dr*x; rx=1/J.array();


// CODE FOR THE CONNECTION OF THE ELEMENTS
int sk,face;
Array2f vn;

vn<<0,1;

SparseMatrix<int> SpFToV(TotalFaces,Nv);
SparseMatrix<int> I(TotalFaces,TotalFaces);
for(int i=0;i<TotalFaces;i++)
I.coeffRef(i,i)=1;

sk=0;
for(int i=0;i<4;i++)
{  for(face=0;face<Nfaces;face++)
{        SpFToV.coeffRef(sk,EToV(i,vn(face)))=1;

sk=sk+1; }}


SpFToV=(SpFToV*SpFToV.transpose());

sk=SpFToV.nonZeros();
sk=sk-TotalFaces;
SpFToV=SpFToV-I;


ArrayXd faces1(sk),faces2(sk),element1(sk),element2(sk),face1(sk),face2(sk);
ArrayXi ind(sk);
int i=0;
for ( k=0;k<SpFToV.outerSize();++k )
    for (SparseMatrix<int>::InnerIterator it(SpFToV,k);it;++it)
    {

        if(it.value()==1)
           { faces1(i)=it.row();
        faces2(i)=it.col();
        i=i+1;
        }
    }

for(i=0;i<sk;i++)

{
element1(i)=ceil((faces1(i)-1)/Nfaces)+1;
element2(i)=ceil((faces2(i)-1)/Nfaces)+1;
   face2(i)=fmod((faces1(i)-1),Nfaces)+1;
   face1(i)=fmod((faces2(i)-1),Nfaces)+1;
   ind(i)=(face1(i)-1)*K+element1(i);
}



EToE.setZero(K*Nfaces,1);
EToF.setZero(K*Nfaces,1);

for(int j=0;j<Nfaces;j++)
for(int i=0;i<K;i++)
{         EToE(i*Nfaces+j,0)=i+1;
     EToF(i*Nfaces+j,0)=j+1;

}
 k=0;
for(int i=0;i<sk;i++)
{
EToE((ind(i)-1),0)=element2(i);
EToF((ind(i)-1),0)=face2(i);
}
EToE.resize(K,Nfaces);
EToF.resize(K,Nfaces);
}





