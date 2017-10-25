#pragma once
#include "common.h"
#include "Interfaces.h"
#include "CommonTIS.h"
#include "linalg.h"


class SCESolver : public IInverseSolver<IFunction> {
public:
	/// Реализация интерфейса
	SolutionPoint solve(const Vector& LowerBound, const Vector& UpperBound, const Vector& X0, IFunction* F, unsigned int MaxCallCount, IInverseSolverEventProcessor* Events);
	void setSearchTerminationCondition(const int LastPointsCount, const double ResidualFuncError);
	TerminationReason getTerminationReason();

	/// Конструктор, деструктор
	SCESolver();
	~SCESolver();

	/// Пользовательские настройки алгоритма
	void setSCEOptions(uint ComplexCount, uint IterCount, uint EvolutionCount, uint MaxFuncEvalCount);

	/// Информация о полученном решении
	SolutionPoint getSolution() const { return currentSolution; } 
	int getFuncEvalCount() const { return nFuncEvals; }
private:
	struct SortData	{
		double Value;
		int Index;

		static bool CompareTo(SortData& A, SortData& B) {
			return A.Value < B.Value;
		}
	};

	typedef vector< vector<double> > Vector2;
	void resize(Vector2& Value, uint N1, uint N2);

	// Генератор случайных чисел
	Random uniRand;

	/// Пользовательские настройки
	uint nComplex;
	uint nIteration;
	uint nEvolution;
	uint maxFuncEvals;

	/// Функция невязки, границы поиска
	IFunction* func;
	uint nDimension;
	Vector x0, lB, uB;
	
	uint nPoints;
	uint nPointsComplex;
	uint nPointsSimplex;
	uint nFuncEvals;
	Vector funcValues;

	/// популяция, комплекс, симплекс ...
	Vector2 complex;
	Vector complexFitness;
	Vector2 population;
	Vector populationFitness;
	Vector2 simplex;
	Vector simplexFitness;
	Vector newPoint;
	double newPointFitness;

	/// вспомогательные переменные
	Vector2 tmpPopulation;
	Vector uBcmx, lBcmx;

	/// текущее решение
	SolutionPoint currentSolution;

	//для условия завершения
	bool doBreak;
	SearchTermination stopper;
	TerminationReason trmnReason;

	/// вывод текущего состояния во вне
	IInverseSolverEventProcessor* events;

	void breakSolve() { doBreak = true; }
	void RunInternal();
	double CalculateCost(Vector const& point);
	void SortPopulation(Vector2& Population, Vector& Population_Fitness );	
	void SortComplex(Vector2& Complex, Vector& Complex_Fitness);
	void CCEUA();
};