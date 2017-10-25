#include "stdafx.h"
#include "common.h"

Vector operator-(const Vector& A, const Vector& B) {
	if (A.size() != B.size()) {
		throw std::exception("Different vector size in Vector::operator-");
	}

	Vector res(A.size());
	for (int i = 0, N = A.size(); i < N; ++i)
		res[i] = A[i] - B[i];

	return res;
}

Vector operator+(const Vector& A, const Vector& B) {
	if (A.size() != B.size()) {
		throw std::exception("Different vector size in Vector::operator-");
	}

	Vector res(A.size());
	for (int i = 0, N = A.size(); i < N; ++i)
		res[i] = A[i] + B[i];

	return res;
}

Vector operator*(double A, const Vector& B) {
	Vector res(B.size());
	for (int i = 0, N = B.size(); i < N; ++i)
		res[i] = A*B[i];

	return res;
}


Vector StringList2Vector(const StringList& Data) {
	Vector line;
	for (int col = 0; col < Data.size(); col++) {
		double v;
		try {
			v = stringTo<double>(Data[col]);
		}
		catch (...) {
			v = std::numeric_limits<double>::quiet_NaN();
		}
		line.push_back(v);
	}
	return line;
}

std::vector<Vector> table2double(const StringTable& Table) {
	std::vector<Vector> res;

	for (int i = 0; i < Table.size(); i++) {
		res.push_back(StringList2Vector(Table[i]));
	}

	return res;
}

bool checkTable4NaN(const std::vector<Vector>& Table) {
	for (int i = 0, N = Table.size(); i < N; ++i) {
		auto v = Table[i];
		for (int j = 0, M = v.size(); j < M; ++j)
				if (_isnan(v[j])) return false;
	}
	return true;
}

// return L2 norm
double fabs(const Vector& X) {
	double sum = 0.0;
	for (int i = 0; i < X.size(); ++i) {
		double Xi = X[i];
		sum += Xi * Xi;
	}
	return sqrt(sum);
}

// return maximum value (module) of vector
double max(const Vector& X) {
	double max_value = 0.0;
	for (int i = 0; i < X.size(); ++i) {
		double Xi = fabs(X[i]);
		max_value = __max(max_value, Xi);
	}
	return max_value;
}

int find(const Vector& X, double Value, double Tolerance) {
	for (int i = 0, N = X.size(); i < N; ++i)
		if (fabs(X[i] - Value) < Tolerance) return i;
	return -1;
}