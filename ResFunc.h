#pragma once

#include "stdafx.h"
#include <fstream>
#include "Interfaces.h"

class ResFunc : IFunction
{
public:
	ResFunc();
	~ResFunc(void);
	Matrix Read(string s, int m, int n);
	double poly(double p, double T);
	double value(const Vector& X); // override
private:
};

