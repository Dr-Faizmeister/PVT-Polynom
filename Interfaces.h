#pragma once

#include "stdafx.h"

/// ��������� �������
class IFunction {
public:
	/// ���������� �������� �������
	virtual double value(const Vector& X) = 0;
};

//typedef boost::shared_ptr<IFunction> Function;

/// ��������� �������� ������� ������
class IGridFunction : public IFunction {
public:
	/// ���������� �������� ��������� ������� ������
	virtual Vector gridValue(const Vector& X) = 0;
	/// �������� ������� ������ ������������, ��� ����� ��������� ���������
	double value(const Vector& X);
};

//typedef boost::shared_ptr<IGridFunction> GridFunction;

/// Interface of 1D function
class IFunction1D {
public:
	/// ���������� �������� �������
	virtual double value(const double x) = 0;
};

/// Reference to the function
//typedef boost::shared_ptr<IFunction1D> Function1D;

/// List of references to the functions
//typedef std::vector<Function1D> Function1DList;

/// ����� ������� ������ ����������� - ������������ ���������� � ���������� ��������
struct SolutionPoint {
	Vector X;
	double Value;

	SolutionPoint(const Vector& _X, double _Value) :
		X(_X), Value(_Value) {}
	SolutionPoint() {}
};

enum TerminationReason {
	TR_MAX_ITERATIONS,
	TR_LAST_POINTS,
	TR_ERROR,
	TR_UNKNOWN,
};

/// ��������� ����������� ������� ��������� ��������
class IInverseSolverEventProcessor {
public:
	/// ���������� ����� ���������� ������� �������� �������
	virtual void OnFunctionCall(const Vector& X, double Value, int FuncEvaluationCout) = 0;
	/// ���������� ����� ���������� ����� ������� ����� �������
	virtual void OnNewSolution(const Vector& X, double Value, int FuncEvaluationCout) = 0;
};

template<class Func> class IInverseSolver {
public:
	/// LowerBound - ������ ������� ������� ��� ������ �������
	/// UpperBound - ������� ������� ������� ��� ������ �������
	/// X0 - ��������� �����������. � ������, ���� ����������� ����� ������� = 0, �������� ������ ��� ������� ��������� ����������� (��� �������������)
	/// F - �������������� �������
	/// MaxCallCount - ����������� ���������� ����� ������� �������. ��� ���������� ����� ����� ������� �������� ������ ������ ���� ��������.
	/// IInverseSolverEventProcessor - ������ ��� ��������� ������� �� ��������� ��������
	/// ����� ���������� ��������� ����� ����������� ��������
	virtual SolutionPoint solve(const Vector& LowerBound, const Vector& UpperBound, const Vector& X0, Func* F, unsigned int MaxCallCount, IInverseSolverEventProcessor* Events) = 0;

	/// LastPointsCount - ���������� ��������� ������ �����, �� ������� ����� �������������� ��������� �������� ������� �������
	/// ResidualFuncError - ��������� ������� �������, ��� ���������� �������� ����� ��������������� ����� ��������
	virtual void setSearchTerminationCondition(const int LastPointsCount, const double ResidualFuncError) = 0;
	virtual TerminationReason getTerminationReason() = 0;
};

