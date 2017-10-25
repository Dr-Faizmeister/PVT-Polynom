#pragma once
#include "stdafx.h"
#include "common.h"

#include <boost/random.hpp>

/// Matrix multiplication
DoubleMatrix MultMatrix(const DoubleMatrix& A, const DoubleMatrix& B);

/// Matrix-vector multiplication
Vector MultMatrixVector(const DoubleMatrix& A, const Vector& B);


/// transpose matrix
DoubleMatrix TransposeMatrix(const DoubleMatrix& A);

/// It solves a system with three-diagonal band matrix by Thomas algorithm in form a*x_{i+1} + b*x_i + c*x_{i-1} = d_i
void Solve1DSLAU(const Vector& a, const Vector& b, const Vector& c, const Vector& d, Vector& x);

/// Gauss linear system solution O(N^3)
void SolveSLAUGauss(const DoubleMatrix& A, const Vector& B, Vector& X);

class TridiagSystem {
public:
	void SetLength(int len);
	void Solve(const Vector& a, const Vector& b, const Vector& c, const Vector& d, Vector& x);
private:
	int N;
	Vector alpha, beta;
};


// Encapsulation of use of boost::random
class Random {
public:
	Random();
	Random(double A, double B);
	Random(double A, double B, uint32_t Seed);
	double operator()() {
		return dist(rng);
	}
private:
	boost::random::mt19937 rng;
	boost::random::uniform_real_distribution<double> dist;
};

typedef DoubleMatrix DoubleArray; // �� �� �����, ��� � DoubleMatrix
typedef std::vector<DoubleArray> MatrixVector;		///< ������ � ������� ������� �����.


/// ������ ���������� ������ �������� NxN
void resize(DoubleArray& A, int N);

/// ������ ������ �������� NxM
void resize(DoubleArray& A, int N, int M);


/// ������ �������� ������ ��� ������-���������������� �������
void ResizeBlockTriDiagSystem(MatrixVector& c, MatrixVector& d, MatrixVector& e, DoubleArray& b, DoubleArray& x, int N, int BlockSize);

/// ������� ������-���������������� ������� �������� ��������� ������� ��������� ��������
void SolveBlockTriDiagSystem(MatrixVector const& c, MatrixVector const& d, MatrixVector const& e, DoubleArray const& b, DoubleArray& x);

// ������� ������-���������������� ������� � �������� ����� 3
// ���������� ������� 3�3 ������ ���� ������������ � ���� ����������� ������� �������� 9 (��������� �� ������� ����� ��)
class C3BlockTridiagSystem {
public:
	void SetNumber(const int Numb);
	void Solve(DoubleArray &a, DoubleArray &c, DoubleArray &b, DoubleArray &d, DoubleArray &x); //a_i x_{i-1} + c_i x_i + b_i x_{i+1} = d_i
private:	
	inline void inverseMatrix(Vector &a, Vector &res);
	inline void calcDenominator(Vector &ci, Vector &ai, Vector &alphai, Vector &res); 
	inline void calcAlpha(Vector &invdeni, Vector &bi, Vector &res);
	inline void calcNumerator(Vector &ai, Vector &di, Vector &betai, Vector &res);  
	inline void calcBeta0(Vector &invdeni, Vector &d0, Vector &res);
	inline void calcBeta(Vector &invdeni, Vector &numi, Vector &res);  
	inline void calcX(Vector &alphai, Vector &xi, Vector &betai, Vector &res);  
	int N;
	DoubleArray alpha, beta;
	Vector den, invden, num;
};

// ������� ������-���������������� ������� � �������� ����� 2
// ���������� ������� 2�2 ������ ���� ������������ � ���� ����������� ������� �������� 9 (��������� �� ������� ����� ��)
class C2BlockTridiagSystem {
public:
	void SetNumber(const int Numb);
	void Solve(DoubleArray &a, DoubleArray &c, DoubleArray &b, DoubleArray &d, DoubleArray &x); //a_i x_{i-1} + c_i x_i + b_i x_{i+1} = d_i
private:	
	inline void inverseMatrix(Vector &a, Vector &res);
	inline void calcDenominator(Vector &ci, Vector &ai, Vector &alphai, Vector &res); 
	inline void calcAlpha(Vector &invdeni, Vector &bi, Vector &res);
	inline void calcNumerator(Vector &ai, Vector &di, Vector &betai, Vector &res);  
	inline void calcBeta0(Vector &invdeni, Vector &d0, Vector &res);
	inline void calcBeta(Vector &invdeni, Vector &numi, Vector &res);  
	inline void calcX(Vector &alphai, Vector &xi, Vector &betai, Vector &res);  
	int N;
	DoubleArray alpha, beta;
	Vector den, invden, num;
};
