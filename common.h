#pragma once
#include "stdafx.h"

// disable warning C4018: '<' : signed/unsigned mismatch
#pragma warning(disable: 4018)

using std::string;
//using boost::make_shared;

typedef unsigned int uint;

typedef std::vector<double> Vector;
typedef std::vector<bool> boolVector;
typedef std::vector<string> StringList;
typedef std::vector< std::vector<double> > Matrix;
typedef std::vector<Vector> DoubleMatrix;

typedef std::vector<string> StringList;
typedef std::vector< std::vector<string> > StringTable;

// return L2 norm
double fabs(const Vector& X);

// return maximum value (module) of vector
double max(const Vector& X);

// find item in Vector
int find(const Vector& X, double Value, double Tolerance = 1e-10);

// True Vertical Depth
struct TVDDepth {
	double tvd;
	explicit TVDDepth(double Value) {
		tvd = Value;
	}
};

typedef std::pair<double, double> DoublePair;
typedef std::vector<DoublePair> DoublePairList;

typedef std::pair<TVDDepth, double> TVDFuncValue;
typedef std::vector<TVDFuncValue> TVDFunc;

/// ===================================================================================
/// Vector tools
/// ===================================================================================

Vector operator-(const Vector& A, const Vector& B);
Vector operator+(const Vector& A, const Vector& B);
Vector operator*(double A, const Vector& B);

/// ===================================================================================
/// String tools
/// ===================================================================================

template<typename T> T stringTo(const std::string& s) {
    return boost::lexical_cast<T>(s);
};

template<typename T> std::string toString(const T& x) {
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << x;
    return s.str();
};

/// It converts the string array to numbers vector. If can not convert the NaN value is inserted
Vector StringList2Vector(const StringList& Data);

/// It converts the string table to numbers table. If can not convert the NaN value is inserted
std::vector<Vector> table2double(const StringTable& Table);

/// Check table of double values for NaN value
bool checkTable4NaN(const std::vector<Vector>& Table);

/// ===================================================================================
/// Data containers and access
/// ===================================================================================

// Осталось после замены своего велосипеда на boost::optional. Семантика привнесена из haskell :)

#define Maybe boost::optional 
#define Nothing() boost::none
#define Just(A) (A)

/// получение последнего элемента вектора
template<typename T> T lastValue(const std::vector<T>& Container) {
	if (Container.size() == 0) {
		throw std::exception("empty container");
	}

	return Container[Container.size() - 1];
}

/// ===================================================================================
/// Normal conditions (атмосферные условия по ГОСТ 2939—63)
/// ===================================================================================
#define NORMAL_CONDITION_PRESSURE 101325.0
#define NORMAL_CONDITION_TEMPERATURE 293.15