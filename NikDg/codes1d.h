#ifndef CODES1D_H
#define CODES1D_H
#include "all.h"

MatrixXd Vandermonde1D(int N,const Ref<const VectorXd>& x);
VectorXd JacobiP(const Ref<const VectorXd>& x,double alpha,double beta,int N);
MatrixXd JacobiGQ(double alpha,double beta,int N);
VectorXd JacobiGL(double alpha,double beta,int N);
VectorXd GradJacobiP(const Ref<const VectorXd>& r,double alpha,double beta,int N);
MatrixXd GradVandermonde1D(int N,const Ref<const VectorXd>& r);
MatrixXd Dmatrix1D(int N,const Ref<const VectorXd>& r,const Ref<const MatrixXd>& V);
MatrixXd GeometricFactors1DJ(const Ref<const MatrixXd>& x,const Ref<const MatrixXd>& Dr);

MatrixXd Connect1DEToE(const Ref<const ArrayXd>& EToV);



#endif // CODES1D_H
