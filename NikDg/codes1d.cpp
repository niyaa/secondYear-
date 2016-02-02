#include "all.h"



VectorXd JacobiP(const Ref<const VectorXd>& x,double alpha,double beta,int N){
      double gamma0, gamma1;
    VectorXd P;
    MatrixXd PL,xp;

    int i,n,j,k;
        n=x.rows();
    PL=MatrixXd::Zero(N+1,n);
    xp=x.transpose();
    gamma0=pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);
    for(k=0;k<n;k++)
        PL(0,k)=1.0/sqrt(gamma0);
    if (N==0)
        return P=PL.row(N);
    gamma1=(alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
    for (k=0;k<n;k++)
       PL(1,k)=((alpha+beta+2)*xp(0,k)/2+(alpha-beta)/2)/sqrt(gamma1);

    if(N==1)
       return P=PL.row(N);
    double aold,anew,bnew,h1;
    aold = 2.0/(2.0+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));
    for (i=0;i<N-1;i++)
    {
        h1=2*(i+1)+alpha+beta;
    anew= 2/(h1+2)*sqrt( (i+2)*(i+2+alpha+beta)*(i+2+alpha)*(i+2+beta)/(h1+1)/(h1+3));
     bnew= -(pow(alpha,2)-pow(beta,2))/h1/(h1+2);

     for (k=0;k<n;k++)
     PL(i+2,k)=1/anew*(-aold*PL(i,k)+(xp(0,k)-bnew)*PL(i+1,k));
aold = anew;
        }
     return P=PL.row(N);
}


MatrixXd Vandermonde1D(int N,const Ref<const VectorXd>& x){
    MatrixXd V1D;
    int j,n;
    n=x.rows();
    V1D=MatrixXd::Zero(n,N+1);
    for (j=0;j<=N;j++)
        V1D.col(j)=JacobiP(x,0,0,j);
   return V1D;
}


MatrixXd JacobiGQ(double alpha,double beta,int N)
{
    MatrixXd xw,J,J1;
    VectorXd h1;
    int i;
xw.setZero(N+1,2);
h1.setZero(N+1);
J1.setZero(N+1,N+1);
if (N==0)
{
   xw(0,0)=-(alpha-beta)/(alpha+beta+2);
   xw(0,1)=2;
   return xw;
}
J.setZero(N+1,N+1);
for(i=0;i<=N;i++)
 {   h1(i)=2*(i)+alpha+beta;
    J1(i,i)=-1/2.0*(pow(alpha,2)-pow(beta,2))/(h1(i)+2.0)/h1(i);
}
for (i=0;i<N;i++)
    J(i,i+1)=2.0/(h1(i)+2)*sqrt((i+1)*((i+1)+alpha+beta)*((i+1)+alpha)*((i+1)+beta)/(h1(i)+1)/(h1(i)+3));
J=J+J1;
if(alpha+beta<2.3e-16)
J(0,0)=0.0;
J1.setZero(N+1,N+1);
J1=J.transpose();
J=J+J1;
//cout<<"the value of J1\n"<<J1<<endl;

SelfAdjointEigenSolver<MatrixXd> eigensolver(J);
if(eigensolver.info()!=Success) abort();

xw.col(0)=eigensolver.eigenvalues();
J1.setZero(N+1,N+1);
J1=eigensolver.eigenvectors();
//cout<<"the eigen vectors are\n"<<J1<<endl;
for(i=0;i<=N;i++)
    xw(i,1)=pow(J1(0,i),2)*pow(2,(alpha,beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);

return xw;
}
VectorXd JacobiGL(double alpha,double beta,int N){
   VectorXd x,y;
    MatrixXd xw;
    int i;
    xw.setZero(N+1,2);
    x.setZero(N+1);
    if(N==1)
    {x(1)=-1.0;
        x(2)=1.0;
     return x;
    }
    xw =JacobiGQ(alpha+1,beta+1,N-2);
    x(0)=-1.0;
    x(N)=1.0;
    y=xw.col(0);
    for(i=1;i<N;i++)
    x(i)=xw(i-1,0);
    return x;
}

VectorXd GradJacobiP(const Ref<const VectorXd>& r,double alpha,double beta,int N){
VectorXd dP;
int n;
n=r.size();
dP.setZero(n);
if(N==0)
    dP(0)=0.0;
else
    dP=sqrt(N*(N+alpha+beta+1))*JacobiP(r,alpha+1,beta+1,N-1);
return dP;

}

MatrixXd GradVandermonde1D(int N,const Ref<const VectorXd>& r){
    // Purpose : Initialize the gradient of the modal basis (i) at (r)  at order N
    MatrixXd DVr;
    int n=r.size();
    DVr.setZero(n,N+1);

    for(int i=0;i<=N;i++)
        DVr.col(i)=GradJacobiP(r,0,0,i);
    return DVr;
}

MatrixXd Dmatrix1D(int N,const Ref<const VectorXd>& r,const Ref<const MatrixXd>& V){
    MatrixXd Vr, Dr;
    Vr=GradVandermonde1D(N,r);
    Dr=Vr*V.inverse();
    //** not sure of the method to solve the linear system **//


    return Dr;
}



MatrixXd Connect1DEToE(const Ref<const ArrayXd>& EToV){

    int Nfaces,K,Nv,TotalFaces,sk,face;

K=EToV.rows();
TotalFaces=Nfaces*K;
Nv=K+1;
Array2f vn;

vn<<0,1;
std::cout<<vn<<"\n";
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
// Total No. of Non Zeros elements, sk
sk=SpFToV.nonZeros();
sk=sk-TotalFaces;
SpFToV=SpFToV-I;


ArrayXd faces1(sk),faces2(sk),element1(sk),element2(sk),face1(sk),face2(sk);
ArrayXi ind(sk);
int i=0;
    for ( int k=0;k<SpFToV.outerSize();++k )
        for (SparseMatrix<int>::InnerIterator it(SpFToV,k);it;++it)
        {
           // cout<<i<<"=\t"<<it.row()<<","<<it.col()<<"\t"<<it.value()<<endl;
            if(it.value()==1)
               { faces1(i)=it.row();
            faces2(i)=it.col();
            i=i+1;
            }
        }

//for(i=0;i<sk;i++)
//    cout<<faces1(i)<<"\t"<<faces2(i)<<endl;



for(i=0;i<sk;i++)

{
    element1(i)=ceil((faces1(i)-1)/Nfaces)+1;
element2(i)=ceil((faces2(i)-1)/Nfaces)+1;
       face2(i)=fmod((faces1(i)-1),Nfaces)+1;
       face1(i)=fmod((faces2(i)-1),Nfaces)+1;
       ind(i)=(face1(i)-1)*K+element1(i);
}

//for(i=0;i<sk;i++)
//std::cout<<ind(i)<<"\n";
ArrayXd EToE,EToF;

EToE.setZero(K*Nfaces);
EToF.setZero(K*Nfaces);

for(int j=0;j<Nfaces;j++)
for(int i=0;i<K;i++)
{         EToE(i*Nfaces+j)=i+1;
         EToF(i*Nfaces+j)=j+1;

}


//cout<<" element2 is"<<element2<<endl;


int k=0;
for(int i=0;i<sk;i++)
{
   EToE(ind(i)-1)=element2(i);
   EToF(ind(i)-1)=face2(i);
}

//cout<<"EToE is"<<EToE<<endl;
//cout<<"EToF is"<<EToF<<endl;



MatrixXd EToEF;
EToEF.setZero(K*Nfaces,2);
EToEF.col(0)=EToE;
EToEF.col(1)=EToF;

//cout<<EToEF<<endl;


    return EToEF;
}

