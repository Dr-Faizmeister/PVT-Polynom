#pragma once
#include "common.h"
#include "Interfaces.h"
#include "CommonTIS.h"
#include "linalg.h"


class SCESolver : public IInverseSolver<IFunction> {
public:
	/// ���������� ����������
	SolutionPoint solve(const Vector& LowerBound, const Vector& UpperBound, const Vector& X0, IFunction* F, unsigned int MaxCallCount, IInverseSolverEventProcessor* Events);
	void setSearchTerminationCondition(const int LastPointsCount, const double ResidualFuncError);
	TerminationReason getTerminationReason();

	/// �����������, ����������
	SCESolver();
	~SCESolver();

	/// ���������������� ��������� ���������
	void setSCEOptions(uint ComplexCount, uint IterCount, uint EvolutionCount, uint MaxFuncEvalCount);

	/// ���������� � ���������� �������
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

	// ��������� ��������� �����
	Random uniRand;

	/// ���������������� ���������
	uint nComplex;
	uint nIteration;
	uint nEvolution;
	uint maxFuncEvals;

	/// ������� �������, ������� ������
	IFunction* func;
	uint nDimension;
	Vector x0, lB, uB;
	
	uint nPoints;
	uint nPointsComplex;
	uint nPointsSimplex;
	uint nFuncEvals;
	Vector funcValues;

	/// ���������, ��������, �������� ...
	Vector2 complex;
	Vector complexFitness;
	Vector2 population;
	Vector populationFitness;
	Vector2 simplex;
	Vector simplexFitness;
	Vector newPoint;
	double newPointFitness;

	/// ��������������� ����������
	Vector2 tmpPopulation;
	Vector uBcmx, lBcmx;

	/// ������� �������
	SolutionPoint currentSolution;

	//��� ������� ����������
	bool doBreak;
	SearchTermination stopper;
	TerminationReason trmnReason;

	/// ����� �������� ��������� �� ���
	IInverseSolverEventProcessor* events;

	void breakSolve() { doBreak = true; }
	void RunInternal();
	double CalculateCost(Vector const& point);
	void SortPopulation(Vector2& Population, Vector& Population_Fitness );	
	void SortComplex(Vector2& Complex, Vector& Complex_Fitness);
	void CCEUA();
};