#include <Eigen/Dense>
#include "mex.h"




using namespace Eigen;




MatrixXcd solver_triang_opt(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(24);
    coeffs[0] = 2*d[0];
    coeffs[1] = d[1];
    coeffs[2] = d[2];
    coeffs[3] = d[3];
    coeffs[4] = d[4];
    coeffs[5] = -d[5];
    coeffs[6] = 2*d[6];
    coeffs[7] = d[7];
    coeffs[8] = d[8];
    coeffs[9] = d[9];
    coeffs[10] = -d[10];
    coeffs[11] = 2*d[11];
    coeffs[12] = d[12];
    coeffs[13] = d[13];
    coeffs[14] = -d[14];
    coeffs[15] = 2*d[15];
    coeffs[16] = d[16];
    coeffs[17] = -d[17];
    coeffs[18] = 2*d[18];
    coeffs[19] = -d[19];
    coeffs[20] = 2;
    coeffs[21] = 1;
    coeffs[22] = -1;
    coeffs[23] = -2;



	// Setup elimination template
	static const int coeffs0_ind[] = { 21,21,21,21,21,22,21,21,20,21,23,22,21,21,4,2,1,3,0,9,7,6,8,1,13,11,7,12,21,2,16,12,8,15,3,18,13,9,16,4,22,21,21,21,21,21,21,21,21,21,21,21,21,22,21,21,22,20,1,0,3,21,20,6,1,8,20,7,2,12,21,20,8,3,15,1,3,0,21,21,23,6,8,1,23,20,9,4,16,7,12,21,2,21,21,21,22,20,2,0,4,7,2,1,1,0,9,7,6,1,4,11,2,2,1,3,0,21,13,9,11,7,7,2,6,8,1,13,11,7,12,2,21,12,2,3,3,1,0,16,12,7,8,3,8,6,1,16,12,11,8,15,3,12,7,2,12,15,8,3,21,13,2,4,4,1,3,0,21,18,13,7,9,9,4,6,8,1,18,13,11,13,9,16,4,7,12,21,2,3,0,21,21,8,21,1,21,12,2,21,21,15,3,22,4,2,21,1,0,3,21,21,9,7,21,6,1,8,10,5,17,13,11,22,21,7,20,2,12,14,5,21,4,2,3,0,1,19,14,10,5,21,9,7,8,1,6,19,14,10,17,5,13,11,12,2,7,22,13,12,16,16,9,4,8,15,3,8,15,3,22,19,14,10,17,5,13,18,22,9,16,4,22,21,23,9,16,4,22,21,16,4,21,21,21,22,23,21 };
	static const int coeffs1_ind[] = { 19,14,17,5,10,14,17,10,5,22,16,12,15,3,8,16,12,22,21,8,3,15,22,14,19,10,17,18,13,16,5,4,9,10,17,18,13,22,23,9,5,4,21,16,17,5,21,22,21,21,19,14,10,22,5,21,17,21,22 };
	static const int C0_ind[] = { 0,57,59,63,118,173,177,230,234,236,291,346,354,359,407,408,413,414,416,465,466,471,472,474,523,524,529,530,531,532,581,582,587,588,590,639,640,645,646,648,649,696,708,767,808,826,869,875,885,944,950,970,1003,1043,1051,1062,1098,1116,1121,1122,1123,1163,1175,1179,1180,1181,1234,1237,1238,1239,1271,1293,1295,1296,1297,1357,1358,1380,1390,1396,1404,1415,1416,1438,1463,1468,1469,1470,1471,1473,1474,1475,1496,1527,1534,1586,1587,1593,1638,1652,1694,1696,1697,1710,1711,1712,1752,1755,1769,1770,1811,1812,1814,1826,1829,1830,1831,1862,1868,1869,1871,1872,1885,1886,1887,1888,1889,1927,1930,1945,1946,1947,1951,1986,1989,2000,2006,2007,2008,2042,2045,2047,2059,2060,2064,2065,2066,2101,2104,2105,2119,2120,2121,2122,2123,2124,2163,2180,2181,2182,2215,2218,2222,2226,2232,2242,2243,2254,2271,2274,2277,2280,2284,2291,2292,2300,2301,2312,2333,2336,2338,2342,2351,2352,2353,2358,2359,2360,2370,2419,2421,2432,2465,2477,2478,2479,2525,2535,2537,2538,2587,2593,2595,2606,2636,2637,2638,2655,2659,2662,2663,2693,2694,2695,2698,2713,2717,2720,2745,2746,2747,2752,2753,2758,2759,2771,2773,2775,2778,2798,2812,2821,2828,2831,2832,2835,2839,2854,2857,2871,2872,2882,2886,2889,2890,2893,2897,2913,2916,2931,2932,2933,2944,2947,2948,2951,2955,2969,2975,2976,2980,2992,2993,2994,2996,2997,3008,3039,3040,3062,3072,3075,3076,3081,3082,3084,3092,3096,3111,3112,3113,3124,3140,3142,3154,3155,3156,3178,3185,3228,3231,3233,3237,3271,3275,3330,3332,3352 } ;
	static const int C1_ind[] = { 44,47,48,51,55,75,92,93,94,95,102,105,106,109,113,142,143,150,152,161,165,168,169,192,196,212,213,218,221,222,224,225,229,255,256,258,259,271,276,277,278,281,282,284,331,333,345,389,391,393,432,433,451,454,455,457,458,513,516 };

	Matrix<double,58,58> C0; C0.setZero();
	Matrix<double,58,9> C1; C1.setZero();
	for (int i = 0; i < 298; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 59; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,58,9> C12 = C0.partialPivLu().solve(C1);




	// Setup action matrix
	// Matrix<double,16, 9> RR;
    MatrixXd RR(16, 9);	
	RR << -C12.bottomRows(7), Matrix<double,9,9>::Identity(9, 9);

	static const int AM_ind[] = { 10,0,1,2,3,4,5,11,6 };
	// Matrix<double, 9, 9> AM;
    MatrixXd AM(9, 9);
	for (int i = 0; i < 9; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 7, 9> sols;
	sols.setZero();

	// Solve eigenvalue problem
	EigenSolver<Matrix<double, 9, 9> > es(AM);
	ArrayXcd D = es.eigenvalues();	
	ArrayXXcd V = es.eigenvectors();

V = (V / V.row(0).array().replicate(9, 1)).eval();


    sols.row(3) = V.row(1).array();
    sols.row(4) = D.transpose().array();
    sols.row(5) = V.row(5).array();
    sols.row(6) = V.row(7).array();




	return sols;
}
// Action =  x5
// Quotient ring basis (V) = 1,x4,x4*x7,x5,x5*x7,x6,x6*x7,x7,x7^2,
// Available monomials (RR*V) = x4*x5,x4*x5*x7,x5^2,x5^2*x7,x5*x6,x5*x6*x7,x5*x7^2,1,x4,x4*x7,x5,x5*x7,x6,x6*x7,x7,x7^2,





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_opt:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_opt:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_opt:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 21 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:triang_opt:incorrectSize", "Input size must be multiple of 21.");
	}
	int n_instances = mxGetNumberOfElements(prhs[0]) / 21;
	double *input = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(7,9*n_instances,mxCOMPLEX);
	double* zr = mxGetPr(plhs[0]);
	double* zi = mxGetPi(plhs[0]);
	for(int k = 0; k < n_instances; k++) {
		const VectorXd data = Map<const VectorXd>(input + k*21, 21);
		MatrixXcd sols = solver_triang_opt(data);
		Index offset = k*sols.size();
		for (Index i = 0; i < sols.size(); i++) {
			zr[i+offset] = sols(i).real();
			zi[i+offset] = sols(i).imag();
		}
	}
}


