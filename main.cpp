// PVT-Polynom.cpp: определяет точку входа для консольного приложения.
//
#pragma once

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "Interfaces.h"
#include "sce.h"
#include "ResFunc.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	extern int n, m;

	const string INPUT_FILE_NAME = "table.txt";
	const string OUTPUT_FILE_NAME = "polynomial.txt";

	ResFunc * rf = new ResFunc();

	rf->Read(INPUT_FILE_NAME, m, n);

	system("pause");
	return 0;
}

