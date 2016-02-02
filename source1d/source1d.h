#ifndef CODES1D_H
#define CODES1D_H
#include"../../usr/include/eigen3/Eigen/Dense"
using namespace Eigen;
MatrixXi meshGen(float xmin,float xmax,int &K,int &Nv,Ref<VectorXd> VX);
MatrixXd Vandermonde1D(int N,const Ref<const VectorXd>& x);
VectorXd JacobiP(const Ref<const VectorXd>& x,double alpha,double beta,int N);
MatrixXd JacobiGQ(double alpha,double beta,int N);
VectorXd JacobiGL(double alpha,double beta,int N);
VectorXd GradJacobiP(const Ref<const VectorXd>& r,double alpha,double beta,int N);
MatrixXd GradVandermonde1D(int N,const Ref<const VectorXd>& r);
MatrixXd Dmatrix1D(int N,const Ref<const VectorXd>& r,const Ref<const MatrixXd>& V);
//VectorXd GeometricFactors1DJrx(const Ref<const VectorXd>& x,const Ref<const MatrixXd>& Dr);
MatrixXd Connect1DEToE(const Ref<const MatrixXi>& EToV);
MatrixXi Normals1D(int Nfp,int Nfaces,int K);


#endif // CODES1D_H

