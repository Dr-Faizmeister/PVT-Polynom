#include "stdafx.h"
#include "Interfaces.h"


double IGridFunction::value(const Vector& X) {
	auto e = this->gridValue(X);
	double sum = 0.0;
	foreach (double ei, e) {
		sum += ei*ei;
	}
	return sum;
}