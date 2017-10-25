#pragma once

#include "stdafx.h"
#include "Interfaces.h"

using namespace std;


double sqr(const double x);

/// Type log data
enum LogType {
	LT_TEMPERATURE,       ///< Temperature
	LT_PRESSURE,          ///< Pressure
};

/// List of log types
typedef std::vector<LogType> LogTypeList;
typedef std::vector<string> StringList;

/// Zone for caclulation of residual function
struct ResidualInterval {
	double begin;               ///> Begin of interval
	double end;                 ///> End of interval
	double weight;              ///> Interval weight for residual function calculation
};

/// List of zones for residual function calculation
typedef std::vector<ResidualInterval> ResidualIntervalList;

/// Table of zones for residual function calculation
typedef std::vector<ResidualIntervalList> ResidualIntervalTable;

/// Table of field data
typedef std::pair<double, double> DoublePair;
typedef std::vector<DoublePair> DoublePairList; 
typedef std::vector<DoublePairList> DoublePairTable;
typedef std::vector<double> Vector;

/// List of integer pairs
typedef std::pair<int, int> IntPair;
typedef std::vector<IntPair> IntPairList;
typedef std::vector<IntPairList> IntPairTable;

/// Piece-linear function
class PieceLinearFunction: public IFunction1D {
public:
	PieceLinearFunction::PieceLinearFunction();
	PieceLinearFunction::PieceLinearFunction(DoublePairList& Points);
	double value(const double x);
	void setPoints(DoublePairList& Points);
private:
	DoublePairList points;
	double x_eps;

	int index(double x);
};

/// End of search condition
class SearchTermination {
public:
	SearchTermination();
	void setConditionValues(const int LastPointsCount, const double ResidualFuncError);
	void addBestFValue(const double bestFValue);
	bool IsTerminated();
private:
	int lastPointsCount;
	double resFuncError;
	vector<double> bestFValues;
};

/// get file name from path
string getFileNameFromFilePath(string path);

/// get path to folder where file is located
string getFolderPathFromFilePath(string path);

/// delete file extension from file path or file name (the symbol "." is not removed)
string deleteFileExtension(string file);

/// exception class for simulator errors
class ESimulatorException {
public:
	ESimulatorException(StringList ErrorMessages): errorMessages(ErrorMessages) {};

	StringList what() {
		return errorMessages;
	};

	string what(int ind) {
		if ( (ind >= 0) && (ind < count()) ) {
			return errorMessages[ind];
		}
		else {
			return "";
		}	
	};

	int count() {
		return errorMessages.size();
	};
private:
	StringList errorMessages;
};