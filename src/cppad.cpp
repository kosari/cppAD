// Sample code for calculation of Jacobian using Automatic Differentiation (AD)
// Author: Sina Nia Kosari
// Download and install CppAD http://www.coin-or.org/CppAD/Doc/download.htm
// Download Eigen from: http://eigen.tuxfamily.org/index.php?title=Main_Page
#include <iostream>
#include <Eigen/Dense>
#include <cppad/cppad.hpp>

using namespace std;
using namespace Eigen;
using CppAD::AD;
using std::vector;

//This function is the target of jacobian: Y = f(t)
// in this example Y0  = t0^2 and Y1 = t0*t1;
// the jacobian for this function is [2*t0 , 0 ; t1 , t0]
template <class T>
T func(T &t1){
	T Y(2);
	Y[0] = t1[0]*t1[0];
	Y[1] = t1[0]*t1[1];
	return Y;
}

//This function performs the jacobian
template <class T>
void jac_func(T (*fp)(T&) , Matrix<double,2,2> &mat){
	vector<AD<double> > X(2),Y(2);
	X[0] = -1e-15;
	X[1] = 0;
	//define the independent variable
	Independent(X);

	Y = fp(X);

	CppAD::ADFun<double> f(X,Y);
	vector<double> jac(2*2),x(2);
	x[0] = -1e-15;
	x[1] = 0;

	//perform the jacobian
	jac = f.Jacobian(x);
	//map to eigen matrix
	Matrix<double,2,2,RowMajor> jac_M = Map<Matrix<double,2,2,RowMajor> >(&jac[0]);
	mat=jac_M;
}

int main()
{
	Matrix<double,2,2> mat= Matrix<double,2,2>::Zero(2,2);
	jac_func<vector<AD<double> > >(func<vector<AD<double> > >,mat);
	cout<<mat;
}

/*
#include <iostream>
#include <Eigen/Dense>
#include <cppad/cppad.hpp>

using namespace std;
using namespace Eigen;
using CppAD::AD;
using std::vector;;

template <class T>
T test_func(T &X){
	T Y(6);
	if (X[0] > 0){
		Y[0] = X[0];
		Y[1] = X[1];
		Y[2] = X[2];
		Y[3] = X[3];
		Y[4] = X[4];
		Y[5] = X[5];
	}
	else{
		Y[0] = -X[0];
		Y[1] = -X[1];
		Y[2] = -X[2];
		Y[3] = -X[3];
		Y[4] = -X[4];
		Y[5] = -X[5];
	}
	return Y;
}

template <class T>
void jacob(T (*func)(T&) , T& point, Matrix<double ,6 ,6> &mat){
	vector<AD<double> > X(6) , Y(6);

	Independent(point);
	Y = func(point);

	CppAD::ADFun<double> f(X,Y);
	vector<double> jac(6*6),x(6);
	for (int i=0 ; i <6 ; i++){
		x[i] = point[i];
	}

	jac = f.Jacobian(x);
	Matrix<double,6,6,RowMajor> jac_M = Map<Matrix<double,6,6,RowMajor> >(&jac[0]);
	mat = jac_M;
}

int main(){
	Matrix<double, 6 , 6> mat = Matrix<double, 6 ,6>::Zero(6,6);
	vector<AD<double> >  point(6);

	for (int i =0 ; i < 6 ; i++){
		point[i] = i;
	}
	jacob<vector<AD<double> > >(test_func<vector<AD<double> > >, point , mat);
	std::cout<<mat<<std::endl;
	return 0;
}
 */
