#include "stdafx.h"
#include "linalg.h"
#include <Windows.h>

DoubleMatrix MultMatrix(const DoubleMatrix& A, const DoubleMatrix& B) {
	if (A.size() == 0 || B.size() == 0) {
		throw std::exception("Empty matrix");
	}
	int K = A[0].size();
	int N = A.size();
	int M = B[0].size();
	if (K != B.size()) {
		throw std::exception("Incompatible matrix sizes");
	}

	DoubleMatrix C;
	resize(C, N, M);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			double sum = 0.0;
			for (int k = 0; k < K; ++k)
				sum += A[i][k]*B[k][j];
			C[i][j] = sum;
		}
	}
	return C;
}

Vector MultMatrixVector(const DoubleMatrix& A, const Vector& B) {
	if (A.size() == 0 || B.size() == 0) {
		throw std::exception("Empty matrix or vector");
	}
	if (A[0].size() != B.size()) {
		throw std::exception("Incompatible matrix-vector sizes");
	}

	Vector C(A.size());
	for (int i = 0; i < A.size(); ++i) {
		double sum = 0.0;
		for (int j = 0; j < B.size(); ++j)
			sum += A[i][j]*B[j];
		C[i] = sum;
	}
	return C;
}

DoubleMatrix TransposeMatrix(const DoubleMatrix& A) {
	int N = A.size();
	int M = A[0].size();

	DoubleMatrix B;
	resize(B, M, N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			B[j][i] = A[i][j];
		}
	}
	return B;
}

void Solve1DSLAU( const Vector& a, const Vector& b, const Vector& c, const Vector& d, Vector& x ) {
	int N = a.size();

	if (N == 1) {
		x[0] = d[0]/b[0];
		return;
	}

	Vector alpha, beta;
	alpha.resize(N); beta.resize(N);

	alpha[1] = -a[0]/b[0];
	beta[1] = d[0]/b[0];
	double b_plus_c_alpha = 0.0;
	for (int i = 1; i < N-1; i++) {
		b_plus_c_alpha = (b[i] + c[i]*alpha[i]);
		alpha[i+1]= -a[i]/b_plus_c_alpha;
		beta[i+1] = (d[i] - c[i]*beta[i])/b_plus_c_alpha;
	}

	x[N-1]=(d[N-1] - c[N-1]*beta[N-1])/(b[N-1]+c[N-1]*alpha[N-1]);
	for (int i = N-1; i > 0; i--) {
		x[i-1] = alpha[i]*x[i] + beta[i];
	}
}

void SolveSLAUGauss(const DoubleMatrix& A, const Vector& B, Vector& X) {
	int N = B.size();
	X.resize(N);
	DoubleMatrix a;
	resize(a, N, N);
	Vector _b(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = A[i][j];
	_b.assign(B.begin(), B.end());

	for (int k = 0; k < (N-1); k++)	{
		for (int i = k+1; i < N; i++) {
			if (fabs(a[k][k]) < 1e-50) 
			{
				throw std::exception("Zero elements on the main diagonal");
			}

			double m = -a[i][k]/a[k][k];
			for (int j = k; j < N; j++)
				a[i][j] += a[k][j]*m;
			_b[i] += _b[k]*m;
		}
	}

	for (int k = N-1; k >= 0; k--) {
		double sum = 0.;
		for (int j = k+1; j < N; j++)
			sum += a[k][j]*X[j];
		X[k] = (_b[k] - sum)/a[k][k];
	}

	Vector resid(N);
	for (int i = 0; i < N; i++)
	{
		resid[i] = B[i];
		for (int j = 0; j < N; j++)
			resid[i] -= A[i][j]*X[j];
	}
}

/////////////////////////////////////////////////////////////////////////

void TridiagSystem::SetLength(int len) {
	N = len;
	alpha.resize(N);
	beta.resize(N);
}

void TridiagSystem::Solve(const Vector& a, const Vector& b, const Vector& c, const Vector& d, Vector& x) {
	if (N == 1) {
		x[0] = d[0]/b[0];
		return;
	}

	double alpha_prev = -a[0]/b[0];
	double beta_prev = d[0]/b[0];
	alpha[1] = alpha_prev;
	beta[1] = beta_prev;
	double b_plus_c_alpha_1 = 0.0;
	for (int i = 1; i < N-1; i++) {
		b_plus_c_alpha_1 = 1.0/(b[i] + c[i]*alpha_prev);
		alpha_prev = -a[i]*b_plus_c_alpha_1;
		beta_prev = (d[i] - c[i]*beta_prev)*b_plus_c_alpha_1;
		alpha[i+1]= alpha_prev;
		beta[i+1] = beta_prev;
	}

	x[N-1]=(d[N-1] - c[N-1]*beta[N-1])/(b[N-1]+c[N-1]*alpha[N-1]);
	for (int i = N-1; i > 0; i--) {
		x[i-1] = alpha[i]*x[i] + beta[i];
	}
}

/////////////////////////////////////////////////////////////////////////
Random::Random() :
	rng(GetTickCount()), dist(0.0, 1.0)
{}

Random::Random(double A, double B) :
	rng(GetTickCount()), dist(A, B)
{}

Random::Random(double A, double B, uint32_t Seed) :
	rng(Seed), dist(A, B)
{}

/////////////////////////////////////////////////////////////////////////

void resize(DoubleArray& A, int N) {
	A.resize(N);
	for (int i = 0; i < N; ++i) {
		A[i].resize(N);
	}
}

void resize(DoubleArray& A, int N, int M) {
	A.resize(N);
	for (int i = 0; i < N; ++i) {
		A[i].resize(M);
	}
}

//////////////////////////////////////////////////////////////////////////
// «десь идет код, необходимый дл€ реализации матричной прогонки
//////////////////////////////////////////////////////////////////////////


//----------------------------------------------------------------------
DoubleArray Minor(DoubleArray const& a, int col, int row) {
	int N = a.size(); 
	int dj, di;
	DoubleArray m(N-1);

	resize(m, N-1);

	for (int i=0; i<N-1; i++)
		for (int j=0; j<N-1; j++)
		{
			if (j<col) dj=0; else dj=1;
			if (i<row) di=0; else di=1;
			m[i][j] = a[i+di][j+dj];
		}

		return m;
}


//----------------------------------------------------------------------
double DetMinor(DoubleArray const& a) {
	int N=a.size();

	switch (N) 
	{
	case  1: 
		return  a[0][0];
	case  2: 
		return  a[0][0]*a[1][1] - a[0][1]*a[1][0];
	case  3: 
		return  a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[1][0]*a[2][1]*a[0][2]
		- a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[1][2]*a[2][1]*a[0][0];
	default:
		{
			double det = 0;
			for (int i=0; i<N; i++)
			{
				if (a[0][i]!=0)
					det += pow(-1.,i)*a[0][i]*DetMinor(Minor(a, i, 0)); 
			}
			return  det;
		}
	} 
}


//----------------------------------------------------------------------
DoubleArray InverseMatrix(DoubleArray const& a) {
	int N = a.size(); 
	DoubleArray x(N);
	DoubleArray ina(N);

	double det = DetMinor(a);

	if (N==1) 
	{
		ina.resize(1);
		ina[0].resize(1); 
		ina[0][0] = 1./det;
	}   
	else
	{
		resize(ina, N);
		resize(x, N);

		if (fabs(det)>0)
		{
			for (int i=0; i<N; i++) {
				Vector& xi = x[i];
				for (int j=0; j<N; j++) {
					xi[j] = pow(-1.,i+j)*DetMinor(Minor(a,j,i))/det;
				}
			}
			for (int i=0; i<N; i++) {
				Vector& ina_i = ina[i];
				for (int j=0; j<N; j++) {
					ina_i[j] = x[j][i];
				}
			}
		}
		else
		{
			for (int i=0; i<N; i++) {
				Vector& ina_i = ina[i];
				const Vector& ai = a[i];
				for (int j=0; j<N; j++)
					ina_i[j] = ai[j];
			}
		}  
	}
	return ina;
}


//----------------------------------------------------------------------
inline Vector ProdMatrixVector(DoubleArray const& a, Vector const& b) {
	int N = a.size(); 
	Vector v;

	v.resize(N);

	for (int i=0; i<N; i++)
	{
		v[i] = 0;
		for (int j=0; j<N; j++)
			v[i] += a[i][j]*b[j];
	}

	return v;
}


//----------------------------------------------------------------------
inline DoubleArray ProdMatrix(DoubleArray const& a, DoubleArray const& b) {
	int N = a.size(); DoubleArray m;

	resize(m, N);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++) {
			m[i][j] = 0;
			for (int k=0; k<N; k++)
				m[i][j] += a[i][k]*b[k][j];
		}

	return m;
}


//----------------------------------------------------------------------
inline DoubleArray SumMatrix(DoubleArray const& a, DoubleArray const& b) {
	int N = a.size(); DoubleArray m;

	resize(m, N);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			m[i][j] = a[i][j] + b[i][j];

	return m;  
}


//----------------------------------------------------------------------
inline Vector SumVector(Vector const& a, Vector const& b) {
	int N = a.size(); Vector v;

	v.resize(N);

	for (int i=0; i<N; i++)
		v[i] = a[i] + b[i];

	return v;  
}


//----------------------------------------------------------------------
inline DoubleArray MinusMatrix(DoubleArray const& a) {
	int N = a.size(); 

	DoubleArray m;
	resize(m, N);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			m[i][j] = - a[i][j];

	return m;  
}


//----------------------------------------------------------------------
inline Vector MinusVector(Vector const& a) {
	int N = a.size(); Vector v;

	v.resize(N);

	for (int i=0; i<N; i++)
		v[i] = - a[i];

	return v;  
}


//----------------------------------------------------------------------
void resizeMatrixVector(MatrixVector& A, int N, int BlockSize) {
	A.resize(N);
	for (int i = 0; i < N; ++i) {
		resize(A[i], BlockSize);
	}
}

void resizeDoubleArray(DoubleArray& A, int N, int BlockSize) {
	A.resize(N);
	for (int i = 0; i < N; ++i) {
		A[i].resize(BlockSize);
	}
}


void ResizeBlockTriDiagSystem(MatrixVector& c, MatrixVector& d, MatrixVector& e, DoubleArray& b, DoubleArray& x, int N, int BlockSize) {
	resizeMatrixVector(c, N, BlockSize);
	resizeMatrixVector(d, N, BlockSize);
	resizeMatrixVector(e, N, BlockSize);
	resizeDoubleArray(b, N, BlockSize);
	resizeDoubleArray(x, N, BlockSize);
}

void SolveBlockTriDiagSystem(MatrixVector const& c, MatrixVector const& d, MatrixVector const& e, DoubleArray const& b, DoubleArray& x) {
	int N = c.size()-1;

	//MatrixVector alpha = c;
	// Ѕыть может имелось в виду?
	MatrixVector alpha;
	alpha.resize(c.size());
	//DoubleArray beta=b;
	DoubleArray beta;
	beta.resize(c.size());

	DoubleArray invD0 = InverseMatrix(d[0]);
	alpha[0] = MinusMatrix(ProdMatrix(invD0, e[0]));
	beta[0]  = ProdMatrixVector(invD0, b[0]);

	for (int i=1; i<N; i++)
	{
		DoubleArray CiAlphai = ProdMatrix(c[i], alpha[i-1]);
		alpha[i] = MinusMatrix(ProdMatrix(InverseMatrix(SumMatrix(d[i], CiAlphai)), e[i]));
		beta[i]  = ProdMatrixVector(InverseMatrix(SumMatrix(d[i], ProdMatrix(c[i], alpha[i-1]))), SumVector(MinusVector(ProdMatrixVector(c[i], beta[i-1])), b[i]));
	}

	x[N] = ProdMatrixVector(InverseMatrix(SumMatrix(d[N], ProdMatrix(c[N], alpha[N-1]))), SumVector(MinusVector(ProdMatrixVector(c[N], beta[N-1])), b[N]));

	for (int i=N-1; i>=0; i--)
		x[i] = SumVector(ProdMatrixVector(alpha[i], x[i+1]), beta[i]);
}


//////////////////////////////////////////////////////////////////////////
// «десь идет код, необходимый дл€ реализации матричной прогонки 
// с заданным размером блока (2 и 3)
//////////////////////////////////////////////////////////////////////////


void C3BlockTridiagSystem::SetNumber(const int Numb) {
	N = Numb;
	alpha.resize(N); beta.resize(N);
	for (int i = 0; i < N; ++i)
	{
		alpha[i].resize(9);
		beta[i].resize(3);
	}
	den.resize(9); invden.resize(9); num.resize(3);
}

//----------------------------------------------------------------------
void C3BlockTridiagSystem::Solve(DoubleArray &a, DoubleArray &c, DoubleArray &b, DoubleArray &d, DoubleArray &x) {
	inverseMatrix(c[0], invden);
	calcAlpha(invden, b[0], alpha[0]);
	calcBeta0(invden, d[0], beta[0]);
	for (int i = 1; i < N; ++i) 	{
		calcDenominator(c[i], a[i], alpha[i-1], den);
		inverseMatrix(den, invden);
		calcAlpha(invden, b[i], alpha[i]); 
		calcNumerator(a[i], d[i], beta[i-1], num);
		calcBeta(invden, num, beta[i]);
	}
	x[N-1] = beta[N-1];
	for (int i = N-2; i >= 0; --i) calcX(alpha[i], x[i+1], beta[i], x[i]);		 
}

//----------------------------------------------------------------------
inline void C3BlockTridiagSystem::inverseMatrix(Vector &a, Vector &res) {
	//с помощью переменных temp просто не считаем два раза повтор€ющиес€ выражени€
	double temp1 = (a[0]*a[8] - a[2]*a[6]);
	double temp2 = (a[5]*a[6] - a[3]*a[8]);
	double temp3 = (a[3]*a[2] - a[5]*a[0]);
	double det = a[4]*temp1 + a[1]*temp2 + a[7]*temp3;
	if (det == 0) assert(false, "Determinant is zero!");
	det = 1./det;
	res[0] = det*(a[4]*a[8] - a[5]*a[7]);
	res[1] = det*(a[2]*a[7] - a[1]*a[8]);
	res[2] = det*(a[1]*a[5] - a[2]*a[4]);
	res[3] = det*temp2;
	res[4] = det*temp1;
	res[5] = det*temp3; 
	res[6] = det*(a[3]*a[7] - a[4]*a[6]);
	res[7] = det*(a[1]*a[6] - a[0]*a[7]);
	res[8] = det*(a[0]*a[4] - a[1]*a[3]);
}

// c[i] + a[i]*alpha[i-1] 
inline void C3BlockTridiagSystem::calcDenominator(Vector &ci, Vector &ai, Vector &alphai, Vector &res) {
	res[0] = ci[0] + ai[0]*alphai[0] + ai[1]*alphai[3] + ai[2]*alphai[6];
	res[1] = ci[1] + ai[0]*alphai[1] + ai[1]*alphai[4] + ai[2]*alphai[7];
	res[2] = ci[2] + ai[0]*alphai[2] + ai[1]*alphai[5] + ai[2]*alphai[8];
	res[3] = ci[3] + ai[3]*alphai[0] + ai[4]*alphai[3] + ai[5]*alphai[6];
	res[4] = ci[4] + ai[3]*alphai[1] + ai[4]*alphai[4] + ai[5]*alphai[7];
	res[5] = ci[5] + ai[3]*alphai[2] + ai[4]*alphai[5] + ai[5]*alphai[8];
	res[6] = ci[6] + ai[6]*alphai[0] + ai[7]*alphai[3] + ai[8]*alphai[6];
	res[7] = ci[7] + ai[6]*alphai[1] + ai[7]*alphai[4] + ai[8]*alphai[7];
	res[8] = ci[8] + ai[6]*alphai[2] + ai[7]*alphai[5] + ai[8]*alphai[8];
}

// - invesred_den[i]*b[i] 
inline void C3BlockTridiagSystem::calcAlpha(Vector &invdeni, Vector &bi, Vector &res) {
	res[0] = - (invdeni[0]*bi[0] + invdeni[1]*bi[3] + invdeni[2]*bi[6]);
	res[1] = - (invdeni[0]*bi[1] + invdeni[1]*bi[4] + invdeni[2]*bi[7]);
	res[2] = - (invdeni[0]*bi[2] + invdeni[1]*bi[5] + invdeni[2]*bi[8]);
	res[3] = - (invdeni[3]*bi[0] + invdeni[4]*bi[3] + invdeni[5]*bi[6]);
	res[4] = - (invdeni[3]*bi[1] + invdeni[4]*bi[4] + invdeni[5]*bi[7]);
	res[5] = - (invdeni[3]*bi[2] + invdeni[4]*bi[5] + invdeni[5]*bi[8]);
	res[6] = - (invdeni[6]*bi[0] + invdeni[7]*bi[3] + invdeni[8]*bi[6]);
	res[7] = - (invdeni[6]*bi[1] + invdeni[7]*bi[4] + invdeni[8]*bi[7]);
	res[8] = - (invdeni[6]*bi[2] + invdeni[7]*bi[5] + invdeni[8]*bi[8]);
}

// d[i] - a[i]*beta[i-1] 
inline void C3BlockTridiagSystem::calcNumerator(Vector &ai, Vector &di, Vector &betai, Vector &res) {
	res[0] = di[0] - (ai[0]*betai[0] + ai[1]*betai[1] + ai[2]*betai[2]);
	res[1] = di[1] - (ai[3]*betai[0] + ai[4]*betai[1] + ai[5]*betai[2]);
	res[2] = di[2] - (ai[6]*betai[0] + ai[7]*betai[1] + ai[8]*betai[2]);
}

// invden*d[0]
inline void C3BlockTridiagSystem::calcBeta0(Vector &invdeni, Vector &d0, Vector &res) {
	res[0] = invdeni[0]*d0[0] + invdeni[1]*d0[1] + invdeni[2]*d0[2];
	res[1] = invdeni[3]*d0[0] + invdeni[4]*d0[1] + invdeni[5]*d0[2];
	res[2] = invdeni[6]*d0[0] + invdeni[7]*d0[1] + invdeni[8]*d0[2];
}

// invden*(- a[i]*beta[i-1] + d[i]) 
inline void C3BlockTridiagSystem::calcBeta(Vector &invdeni, Vector &numi, Vector &res) {
	res[0] = invdeni[0]*numi[0] + invdeni[1]*numi[1] + invdeni[2]*numi[2];
	res[1] = invdeni[3]*numi[0] + invdeni[4]*numi[1] + invdeni[5]*numi[2];
	res[2] = invdeni[6]*numi[0] + invdeni[7]*numi[1] + invdeni[8]*numi[2];
}

// alpha[i]*x[i+1] + beta[i] 
inline void C3BlockTridiagSystem::calcX(Vector &alphai, Vector &xi, Vector &betai, Vector &res) {
	res[0] = alphai[0]*xi[0] + alphai[1]*xi[1] + alphai[2]*xi[2] + betai[0];
	res[1] = alphai[3]*xi[0] + alphai[4]*xi[1] + alphai[5]*xi[2] + betai[1];
	res[2] = alphai[6]*xi[0] + alphai[7]*xi[1] + alphai[8]*xi[2] + betai[2];
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void C2BlockTridiagSystem::SetNumber(const int Numb) {
	N = Numb;
	alpha.resize(N); beta.resize(N);
	for (int i = 0; i < N; ++i) {
		alpha[i].resize(4);
		beta[i].resize(2);
	}
	den.resize(4); invden.resize(4); num.resize(2);
}

//----------------------------------------------------------------------
void C2BlockTridiagSystem::Solve(DoubleArray &a, DoubleArray &c, DoubleArray &b, DoubleArray &d, DoubleArray &x) {
	inverseMatrix(c[0], invden);
	calcAlpha(invden, b[0], alpha[0]);
	calcBeta0(invden, d[0], beta[0]);
	for (int i = 1; i < N; ++i)
	{
		calcDenominator(c[i], a[i], alpha[i-1], den);
		inverseMatrix(den, invden);
		calcAlpha(invden, b[i], alpha[i]); 
		calcNumerator(a[i], d[i], beta[i-1], num);
		calcBeta(invden, num, beta[i]);
	}
	x[N-1] = beta[N-1];
	for (int i = N-2; i >= 0; --i) calcX(alpha[i], x[i+1], beta[i], x[i]);		 
}

//----------------------------------------------------------------------
inline void C2BlockTridiagSystem::inverseMatrix(Vector &a, Vector &res) {
	double det = a[0]*a[3] - a[1]*a[2];
	if (det == 0) assert(false, "Determinant is zero!");
	det = 1./det;
	res[0] = det*a[3];
	res[1] = - det*a[1];
	res[2] = - det*a[2];
	res[3] = det*a[0];
}

// c[i] + a[i]*alpha[i-1] 
inline void C2BlockTridiagSystem::calcDenominator(Vector &ci, Vector &ai, Vector &alphai, Vector &res) {
	res[0] = ci[0] + ai[0]*alphai[0] + ai[1]*alphai[2];
	res[1] = ci[1] + ai[0]*alphai[1] + ai[1]*alphai[3];
	res[2] = ci[2] + ai[2]*alphai[0] + ai[3]*alphai[2];
	res[3] = ci[3] + ai[2]*alphai[1] + ai[3]*alphai[3];
}

// - invesred_den[i]*b[i] 
inline void C2BlockTridiagSystem::calcAlpha(Vector &invdeni, Vector &bi, Vector &res) {
	res[0] = - (invdeni[0]*bi[0] + invdeni[1]*bi[2]);
	res[1] = - (invdeni[0]*bi[1] + invdeni[1]*bi[3]);
	res[2] = - (invdeni[2]*bi[0] + invdeni[3]*bi[2]);
	res[3] = - (invdeni[2]*bi[1] + invdeni[3]*bi[3]);
}

// d[i] - a[i]*beta[i-1] 
inline void C2BlockTridiagSystem::calcNumerator(Vector &ai, Vector &di, Vector &betai, Vector &res) {
	res[0] = di[0] - (ai[0]*betai[0] + ai[1]*betai[1]);
	res[1] = di[1] - (ai[2]*betai[0] + ai[3]*betai[1]);
}

// invden*d[0]
inline void C2BlockTridiagSystem::calcBeta0(Vector &invdeni, Vector &d0, Vector &res) {
	res[0] = invdeni[0]*d0[0] + invdeni[1]*d0[1];
	res[1] = invdeni[2]*d0[0] + invdeni[3]*d0[1];
}

// invden*(- a[i]*beta[i-1] + d[i]) 
inline void C2BlockTridiagSystem::calcBeta(Vector &invdeni, Vector &numi, Vector &res) {
	res[0] = invdeni[0]*numi[0] + invdeni[1]*numi[1];
	res[1] = invdeni[2]*numi[0] + invdeni[3]*numi[1];
}

// alpha[i]*x[i+1] + beta[i] 
inline void C2BlockTridiagSystem::calcX(Vector &alphai, Vector &xi, Vector &betai, Vector &res) {
	res[0] = alphai[0]*xi[0] + alphai[1]*xi[1] + betai[0];
	res[1] = alphai[2]*xi[0] + alphai[3]*xi[1] + betai[1];
}
