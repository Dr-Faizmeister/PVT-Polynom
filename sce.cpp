#include "stdafx.h"
#include <limits>
#include <math.h>
#include <algorithm>

#include "sce.h"
#include "linalg.h"
#include "MessagesLogs.h"

SCESolver::SCESolver() : doBreak(false), trmnReason(TR_UNKNOWN), uniRand(0.0, 1.0) { 
	currentSolution.Value = 100e+100; 

	nComplex = 5;             //number of complexes  (in article this parameter is called "p")
	nIteration = 20;          //number of internal iterations (in article this parameter is called "beta")
	nEvolution = 1;           //number of evolutions of each complex (in article this parameter is called "alpha")
	maxFuncEvals = 25000;     //default number of maximum function evaluations (this parameter could be set in TIS input file)
}

SCESolver::~SCESolver() {
}

SolutionPoint SCESolver::solve(const Vector& LowerBound, const Vector& UpperBound, const Vector& X0, 
							   IFunction* Func, unsigned int MaxCallCount, IInverseSolverEventProcessor* Events) {
	events = Events;
	maxFuncEvals = MaxCallCount;

	func = Func;
	x0 = X0;
	lB = LowerBound;
	uB = UpperBound;

	RunInternal();

	return currentSolution;
}

void SCESolver::setSearchTerminationCondition(const int LastPointsCount, const double ResidualFuncError) {
	stopper.setConditionValues(LastPointsCount, ResidualFuncError);
}

TerminationReason SCESolver::getTerminationReason() {
	return trmnReason;
}

void SCESolver::RunInternal() {
	funcValues.clear();

	// проверка входных параметров
	if (lB.size() != x0.size()) throw std::exception("LB and X0 have incompatible dimensions in SCE");
	if (uB.size() != x0.size()) throw std::exception("UB and X0 have incompatible dimensions in SCE");

	// Число переменных, подлежащих оптимизации
	nDimension = x0.size();

	// Step 0 begin 
	nPointsComplex = 2 * nDimension + 1;
	nPointsSimplex = nDimension + 1;
	nPoints = nComplex * nPointsComplex;
	// Step 0 end

	// Step 1 begin
	// Инициализация популяции. Это - двумерный массив, первый индекс которого обозначает номер точки, второй - ее координата
	resize(population, nPoints, nDimension);
	resize(tmpPopulation, nPoints, nDimension);
	// Инициализация массива для хранения комплекса
	resize(complex, nPointsComplex, nDimension);
	complexFitness.resize(nPointsComplex);
	// Инициализация массива для хранения симлплекса
	resize(simplex, nPointsSimplex, nDimension);
	simplexFitness.resize(nPointsSimplex);
	//
	lBcmx.resize(nDimension);
	uBcmx.resize(nDimension);
	newPoint.resize(nDimension);

	// добавляем в популяцию пользовательское начальное приближение
	population[0] = x0;
	
	// остальные начальные точки популяции рассчитываем исходя из равномерного распределения по области поиска
	for (int i = 1; i < nPoints; i++) {
		for (int j = 0; j < nDimension; j++) {
			population[i][j] = lB[j] + uniRand() * (uB[j] - lB[j]);
		}
	}

	uint nIterations(0);
	nFuncEvals = 0;

	// подсчитываем стоимость каждой точки в начальной популяции
	populationFitness.resize(nPoints);
	for (int i = 0; i < nPoints; i++) {
		if (doBreak) {
			break;
		}
		populationFitness[i] = CalculateCost(population[i]);
	}
	// Step 1 end

	// Step 2 begin
	// Сортируем популяцию в порядке возрастания стоимости точки
	SortPopulation(population, populationFitness);
	// Step 2 end

	// Для каждой итерации
	do {
		if (doBreak) {
			break;
		}

		nIterations++;

		// Исходная популяция разбивается на комплексы
		// для каждого комплекса...
		for (int j = 0; j < nComplex; j++) {
			if (doBreak) break;

			// Step 3 begin
			// Конструируем j-й комплекс из популяции
			for (int k1 = 0; k1 < nPointsComplex; k1++) {
				int k2 = k1 * nComplex + j;
				complex[k1] = population[k2];
				complexFitness[k1] = populationFitness[k2];
			}
			// Step 3 end

			// Step 4 begin
			// С каждым комплексом что-то делаем согласно competitive complex evolution (CCE) algorithm as described in Duan et al. (1992)
			for (int k = 0; k < nIteration; k++) {
				if (doBreak) {
					break;
				}

				// Substeps 1 and 2 (inside Step 4) begin
				vector<int> location;

				for (int l = 0; l < nPointsSimplex; l++) {
					if (doBreak) {
						break;
					}

					int Dummy = std::numeric_limits<int>::min();
					bool Already_Parent = false;
					while (!Already_Parent)
					{
						if (doBreak) break;
						//Dummy = (int)(floor(nPoints_Complex + 0.5 - sqrt(pow(nPoints_Complex + 0.5, 2) - nPoints_Complex * (nPoints_Complex + 1) * norm_rand())));
						Dummy = floor( nPointsComplex + 0.5 - sqrt (pow(nPointsComplex + 0.5, 2) - nPointsComplex * (nPointsComplex + 1) * uniRand() ) );
						//в выражении выше можно получить Dummy = nPoints_Complex (хоть и очень маловероятно) и вылететь за пределы комплекса
						if (Dummy >= nPointsComplex)
							Dummy = nPointsComplex - 1;
						Already_Parent = (std::find(location.begin(), location.end(), Dummy) == location.end());
					}
					location.push_back(Dummy);
				}
				// Substeps 1 and 2 (inside Step 4) end		
				
				// Substep 3 (inside Step 4) begin
				std::sort(location.begin(), location.end());

				// Конструируем симплекс
				for (int l = 0; l < nPointsSimplex; l++) {
					if (doBreak) {
						break;
					}

					simplex[l] = complex[location[l]];
					simplexFitness[l] = complexFitness[location[l]];
				}
						
				// Генерируем новые точки для симплекса
				newPointFitness = std::numeric_limits<double>::quiet_NaN();

				//вычисляем минимальный гиперкуб, содержащий комплекс
				lBcmx = complex[0];
				uBcmx = complex[0];

				for (int i = 1; i < nPointsComplex; i++) {
					for (int id = 0; id < nDimension; id++) {
						if (lBcmx[id] > complex[i][id]) {
							lBcmx[id] = complex[i][id];
						}
						if (uBcmx[id] < complex[i][id]) {
							uBcmx[id] = complex[i][id];
						}
					}
				}

				for (int i = 0; i < nEvolution; i++) {
					CCEUA();

					// Замена худшей точки на некий результат CCEUA?
					simplex[nPointsSimplex - 1] = newPoint;
					simplexFitness[nPointsSimplex - 1] = newPointFitness;
				}
				// Substep 3 (inside Step 4) end

				// Substed 4 (inside Step 4) begin
				// replace the simplex into the complex
				for (int l = 0; l < nPointsSimplex; l++) {
					complex[location[l]] = simplex[l];
					complexFitness[location[l]] = simplexFitness[l];
				}

				SortComplex(complex, complexFitness);
				// Substep 4 (inside Step 4) end

				// Step 4 end
			}

			// Step 5 begin
			// replace the complex back into the population
			for (uint k1 = 0; k1 < nPointsComplex; k1++) {
				int k2 = k1 * nComplex + j;

				population[k2] = complex[k1];
				populationFitness[k2] = complexFitness[k1];
			}
		}

		SortPopulation(population, populationFitness);
		// Step 5 end
	}
	while (!doBreak);

}

void SCESolver::resize( Vector2& Value, uint N1, uint N2 ) {
	Value.resize(N1);
	foreach(Vector& v, Value) v.resize(N2);
}

void SCESolver::SortComplex( Vector2& Complex, Vector& ComplexFitness ) {
	vector<SortData> list;
	for (uint i = 0; i < ComplexFitness.size(); i++) {
		SortData data;
		data.Value = ComplexFitness[i];
		data.Index = i;

		list.push_back(data);
	}
	std::sort(list.begin(), list.end(), SortData::CompareTo);

	Vector2 tmpComplex = Complex;

	for (uint i = 0; i < ComplexFitness.size(); i++) {
		ComplexFitness[i] = list[i].Value;
		Complex[i] = tmpComplex[list[i].Index];
	}
}

void SCESolver::SortPopulation( Vector2& Population, Vector& PopulationFitness ) {
	vector<SortData> list;
	for (uint i = 0; i < nPoints; i++) {
		SortData data;
		data.Value = PopulationFitness[i];
		data.Index = i;

		list.push_back(data);
	}
	std::sort(list.begin(), list.end(), SortData::CompareTo);

	tmpPopulation = Population;

	for (uint i = 0; i < nPoints; i++) {
		PopulationFitness[i] = list[i].Value;
		Population[i] = tmpPopulation[list[i].Index];
	}
}

void SCESolver::CCEUA() {
	SortComplex(simplex, simplexFitness);
	double oldPointFitness = simplexFitness[nPointsSimplex - 1];

	// Вычисляем центроид
	Vector ce(nDimension);
	for (int i = 0; i < nDimension; i++) {
		ce[i] = 0;
		for (int j = 0; j < nPointsSimplex-1; j++) {
			ce[i] += simplex[j][i];
		}
		ce[i] /= (nPointsSimplex - 1.0);
	}

	// Новая точка (отражение)
	for (int i = 0; i < nDimension; i++) {
		newPoint[i] = 2.0 * ce[i] - simplex[nPointsSimplex - 1][i]; 
	}

	// Проверка на выход за границы диапазона
	bool isOutOfBounds = false;
	for (int i = 0; i < nDimension; i++) {
		if (newPoint[i] < lB[i] || newPoint[i] > uB[i]) {
			isOutOfBounds = true;
			break;
		}
	}

	// мутация
	if (isOutOfBounds) {
		for (int i = 0; i < nDimension; i++) {
			newPoint[i] = lBcmx[i] + uniRand() * (uBcmx[i] - lBcmx[i]);
		}
	}

	newPointFitness = CalculateCost(newPoint);

	if (newPointFitness > oldPointFitness) {
		// Новая точка (сжатие)
		for (int i = 0; i < nDimension; i++) {
			newPoint[i] = 0.5 * (ce[i] + simplex[nPointsSimplex - 1][i]);
		}

		newPointFitness = CalculateCost(newPoint);

		if (newPointFitness > oldPointFitness) {
			for (int i = 0; i < nDimension; i++) {
				newPoint[i] = lBcmx[i] + uniRand()  * (uBcmx[i] - lBcmx[i]);
			}
			newPointFitness = CalculateCost(newPoint);
		}
	}
}

double SCESolver::CalculateCost( Vector const& point) {
	if (nFuncEvals > maxFuncEvals) {
		trmnReason = TR_MAX_ITERATIONS;
		breakSolve();
	}

	// Здесь присутствуют некие магические числа. Суть их состоит в том, чтобы при вылете точки за границы куба, в котором проводится минимизация,
	// значение функции получилось бы очень большим числом.
	for (uint i = 0; i < nDimension; i++)	{
		if (point[i] < lB[i]) {
			WMSG << LOG_DEBUG_INFO << "Out of parameters range: param[" << (int)i << "]=" << point[i] << "\n";
			return 1e+12 + (lB[i] - point[i]) * 1e+6;
		}

		if (point[i] > uB[i]) {
			WMSG << LOG_DEBUG_INFO << "Out of parameters range: param[" << (int)i << "]=" << point[i] << "\n";
			return 1e+12 + (point[i] - uB[i]) * 1e+6;
		}
	}

	nFuncEvals++;
	double fValue = func->value(point);

	//проверка, корректно ли решение
	if ( (fValue != fValue) || ( (fValue - fValue) != (fValue - fValue) ) ) {
		trmnReason = TR_ERROR;
		doBreak = true;
	}	

	funcValues.push_back(fValue);

	events->OnFunctionCall(point, fValue, nFuncEvals);
	if (fValue < currentSolution.Value) {
		currentSolution.Value = fValue;
		currentSolution.X = point;

		events->OnNewSolution(point, fValue, nFuncEvals);

		stopper.addBestFValue(fValue);
		if (stopper.IsTerminated()) {
			trmnReason = TR_LAST_POINTS;
			breakSolve();
		}
	}

	return fValue;
}

