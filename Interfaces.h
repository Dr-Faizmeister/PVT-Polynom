#pragma once

#include "stdafx.h"

/// Интерфейс функции
class IFunction {
public:
	/// Возвращает значение функции
	virtual double value(const Vector& X) = 0;
};

//typedef boost::shared_ptr<IFunction> Function;

/// Интерфейс сеточной функции ошибки
class IGridFunction : public IFunction {
public:
	/// Возвращает значение компонент функции ошибки
	virtual Vector gridValue(const Vector& X) = 0;
	/// Значение функции ошибки определяется, как сумма квадратов компонент
	double value(const Vector& X);
};

//typedef boost::shared_ptr<IGridFunction> GridFunction;

/// Interface of 1D function
class IFunction1D {
public:
	/// Возвращает значение функции
	virtual double value(const double x) = 0;
};

/// Reference to the function
//typedef boost::shared_ptr<IFunction1D> Function1D;

/// List of references to the functions
//typedef std::vector<Function1D> Function1DList;

/// Точка решения задачи минимизации - совокупность параметров и найденного значения
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

/// Интерфейс обработчика событий обратного решателя
class IInverseSolverEventProcessor {
public:
	/// Вызывается после очередного расчета значения функции
	virtual void OnFunctionCall(const Vector& X, double Value, int FuncEvaluationCout) = 0;
	/// Вызывается после нахождения более хорошей точки решения
	virtual void OnNewSolution(const Vector& X, double Value, int FuncEvaluationCout) = 0;
};

template<class Func> class IInverseSolver {
public:
	/// LowerBound - нижняя граница области для поиска решения
	/// UpperBound - верхняя граница области для поиска решения
	/// X0 - начальное приближение. В случае, если размерность этого вектора = 0, алгоритм должен сам выбрать начальное приближение (при необходимости)
	/// F - минимизируемая функция
	/// MaxCallCount - максимально допустимое число вызовов функции. При превышении этого числа решение обратной задачи должно быть прервано.
	/// IInverseSolverEventProcessor - объект для обработки событий от обратного решателя
	/// Метод возвращает найденную точку глобального минимума
	virtual SolutionPoint solve(const Vector& LowerBound, const Vector& UpperBound, const Vector& X0, Func* F, unsigned int MaxCallCount, IInverseSolverEventProcessor* Events) = 0;

	/// LastPointsCount - количество последних лучших точек, по которым будет рассчитываться изменение значения функции невязки
	/// ResidualFuncError - изменение функции невязки, при достижении которого будет останавливаться поиск минимума
	virtual void setSearchTerminationCondition(const int LastPointsCount, const double ResidualFuncError) = 0;
	virtual TerminationReason getTerminationReason() = 0;
};

