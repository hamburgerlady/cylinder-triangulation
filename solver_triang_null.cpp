#include <Eigen/Dense>
#include "mex.h"




using namespace Eigen;




MatrixXcd solver_triang_null(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(12);
    coeffs[0] = -d[5]*d[6] + 2*d[2]*d[8] + d[5]*d[9] - 2*d[4]*d[10] - d[0]*d[11] + d[3]*d[11];
    coeffs[1] = -d[5]*d[12] + 2*d[2]*d[14] + d[5]*d[15] - 2*d[4]*d[16] - d[0]*d[17] + d[3]*d[17];
    coeffs[2] = -d[11]*d[12] + 2*d[8]*d[14] + d[11]*d[15] - 2*d[10]*d[16] - d[6]*d[17] + d[9]*d[17];
    coeffs[3] = std::pow(d[2],2) - std::pow(d[4],2) - d[0]*d[5] + d[3]*d[5];
    coeffs[4] = -d[5]*d[7] + d[4]*d[8] + d[2]*d[10] - d[1]*d[11];
    coeffs[5] = -d[5]*d[13] + d[4]*d[14] + d[2]*d[16] - d[1]*d[17];
    coeffs[6] = std::pow(d[8],2) - std::pow(d[10],2) - d[6]*d[11] + d[9]*d[11];
    coeffs[7] = -d[11]*d[13] + d[10]*d[14] + d[8]*d[16] - d[7]*d[17];
    coeffs[8] = std::pow(d[14],2) - std::pow(d[16],2) - d[12]*d[17] + d[15]*d[17];
    coeffs[9] = d[2]*d[4] - d[1]*d[5];
    coeffs[10] = d[8]*d[10] - d[7]*d[11];
    coeffs[11] = d[14]*d[16] - d[13]*d[17];



	// Setup elimination template
	static const int coeffs0_ind[] = { 9,3,4,9,3,0,10,4,0,6,5,3,9,1,7,5,1,0,4,2,10,6 };
	static const int coeffs1_ind[] = { 8,11,11,1,5,8,11,8,2,7,7,2,6,10 };
	static const int C0_ind[] = { 0,5,6,7,8,11,12,13,14,17,18,21,22,23,24,25,26,27,28,29,31,32 } ;
	static const int C1_ind[] = { 3,4,6,9,10,11,13,14,15,16,19,20,21,22 };

	Matrix<double,6,6> C0; C0.setZero();
	Matrix<double,6,4> C1; C1.setZero();
	for (int i = 0; i < 22; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 14; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,6,4> C12 = C0.partialPivLu().solve(C1);



	// Setup action matrix
	// Matrix<double,6, 4> RR;
    MatrixXd RR(6, 4);	
	RR << -C12.bottomRows(2), Matrix<double,4,4>::Identity(4, 4);

	static const int AM_ind[] = { 4,0,5,1 };
	// Matrix<double, 4, 4> AM;
    MatrixXd AM(4, 4);
	for (int i = 0; i < 4; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 2, 4> sols;
	sols.setZero();

	// Solve eigenvalue problem
	EigenSolver<Matrix<double, 4, 4> > es(AM);
	ArrayXcd D = es.eigenvalues();	
	ArrayXXcd V = es.eigenvectors();

V = (V / V.row(0).array().replicate(4, 1)).eval();


    sols.row(0) = V.row(1).array();
    sols.row(1) = D.transpose().array();









	return sols;
}
// Action =  y
// Quotient ring basis (V) = 1,x,y,y^2,
// Available monomials (RR*V) = x*y,y^3,1,x,y,y^2,





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_null:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_null:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_null:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 18 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_null:incorrectSize", "Input size must be multiple of 18.");
	}
	int n_instances = mxGetNumberOfElements(prhs[0]) / 18;
	double *input = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(2,4*n_instances,mxCOMPLEX);
	double* zr = mxGetPr(plhs[0]);
	double* zi = mxGetPi(plhs[0]);
	for(int k = 0; k < n_instances; k++) {
		const VectorXd data = Map<const VectorXd>(input + k*18, 18);
		MatrixXcd sols = solver_triang_null(data);
		Index offset = k*sols.size();
		for (Index i = 0; i < sols.size(); i++) {
			zr[i+offset] = sols(i).real();
			zi[i+offset] = sols(i).imag();
		}
	}
}


