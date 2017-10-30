// PVT-Polynom.cpp: определяет точку входа для консольного приложения.
//
#pragma once

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "Interfaces.h"
#include "sce.h"
#include "ResFunc.h"
#include "EventProc.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	const int COLUMNS = 11;
	const int ROWS = 7;
	const string INPUT_FILE_NAME = "table.txt";
	const string OUTPUT_FILE_NAME = "polynomial.txt";
	Vector LOWER_BOUND(6);
	Vector UPPER_BOUND(6);
	Vector INITIAL_VALUE(6);
	unsigned int MAX_CALL_COUNT = 1000;
	const int LAST_POINTS_COUNT = 10;
	const int RESIDUAL_FUNC_ERROR = 1e-10;
	
	for (double lb : LOWER_BOUND) { lb = -1e10; }

	for (double ub : UPPER_BOUND) { ub = 1e10; }

	for (double iv : INITIAL_VALUE) { iv = 0; }

	ResFunc * rf = new ResFunc();

	rf->Read(INPUT_FILE_NAME, ROWS, COLUMNS);

	SCESolver * is = new SCESolver();

	is->setSearchTerminationCondition(LAST_POINTS_COUNT, RESIDUAL_FUNC_ERROR);
	
	EventProc * events = new EventProc();

	is->solve(LOWER_BOUND, UPPER_BOUND, INITIAL_VALUE, * rf, MAX_CALL_COUNT, * events);

	system("pause");
	return 0;
}

