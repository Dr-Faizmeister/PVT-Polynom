#include "stdafx.h"
#include "CommonTIS.h"


double sqr(const double x) {
	return x*x;
};

//----------------------------------------------------------------
//-----------------PieceLinearFunction---------------------------
//---------------------------------------------------------------

PieceLinearFunction::PieceLinearFunction() {
	x_eps = 1e-10;
};

PieceLinearFunction::PieceLinearFunction(DoublePairList& Points): points(Points) {
	x_eps = 1e-10;
}

int PieceLinearFunction::index(double x) {
	//пытаемся угадать индекс нужной точки, исходя из предположения о том, что шаг постоянный (dX)
	int count = points.size(); 
	int apprInd = int((count - 1)*(x - points.front().first)/(points.back().first - points.front().first));
	if (apprInd == (count - 1)) {
		--apprInd;
	}

	if ((x >= points[apprInd].first - x_eps) && (x <= points[apprInd+1].first + x_eps)) {
		return apprInd;
	}

    if (x > points[apprInd + 1].first) {
		for (int i = apprInd + 1; i < count - 1; i++) {
			if ((x >= points[i].first - x_eps) && (x <= points[i+1].first + x_eps)) {
				return i;
			}
		}
	}

	if (x < points[apprInd].first) {
		for (int i = apprInd; i > 0; i--) {
			if ((x >= points[i-1].first - x_eps) && (x <= points[i].first + x_eps)) {
				return i - 1;
			}
		}
	}

	throw std::exception(("Argument " + toString(x) + " is out of function range: " + toString(points.front().first) + ", " + toString(points.back().first) + " (index)").c_str());
}

double PieceLinearFunction::value(const double x) {
	if ((x < points.front().first - 0.1*x_eps) || (x > points.back().first + 0.1*x_eps)) {
		throw std::exception(("Argument " + toString(x) + " is out of function range: " + toString(points.front().first) + ", " + toString(points.back().first) + " (value)").c_str());
	}

	int ind = index(x);

	if (fabs(points[ind+1].first - points[ind].first) < x_eps) {
		return 0.5*(points[ind+1].second + points[ind].second);
	}

	return points[ind].second + (points[ind+1].second - points[ind].second) * (x - points[ind].first)/(points[ind+1].first - points[ind].first);
}

void PieceLinearFunction::setPoints(DoublePairList& Points) {
	points = Points;
}


//----------------------------------------------------------------
//-----------------End of search condition-----------------------
//---------------------------------------------------------------

SearchTermination::SearchTermination() {
	lastPointsCount = 10; 
	resFuncError = 1e-15;
}

void SearchTermination::setConditionValues(const int LastPointsCount, const double ResidualFuncError) {
	lastPointsCount = LastPointsCount; 
	resFuncError = ResidualFuncError;
}

void SearchTermination::addBestFValue(const double bestFValue) {
	if (bestFValues.size() == 0) {
		bestFValues.push_back(bestFValue);
	}
	else {
		double lastF = bestFValues[bestFValues.size() - 1];
		if (bestFValue < lastF) {
			bestFValues.push_back(bestFValue);
		}
	}
}

bool SearchTermination::IsTerminated() {
	int count = bestFValues.size();

	if ((count < lastPointsCount) || (lastPointsCount <= 1)) {
		return false;
	}

	int bgnInd = count - lastPointsCount;
	return fabs(bestFValues[bgnInd] - bestFValues[count-1]) <= resFuncError;
}

//----------------------------------------------------------------
//-----------------Some functions for work with file path--------
//--------------------------------------------------------------

string getFileNameFromFilePath(string path) {
	//получаем перевернутое имя 
	string fileNameRev = "";
	while ( (!path.empty()) && (path.back() != '/') && (path.back() != '\\') ) {
		fileNameRev.push_back(path.back());	
		path.pop_back();
	}

	//разворачиваем имя
	string fileName = "";
	for (int i = fileNameRev.size() - 1; i >= 0; i--) {
		fileName.push_back(fileNameRev[i]);
	}

	return fileName;
}

string getFolderPathFromFilePath(string path) {
	while ((!path.empty()) && (path.back() != '/') && (path.back() != '\\')) {
		path.pop_back();
	}

	return path;
}

string deleteFileExtension(string file) {
	while ( (!file.empty()) && (file.back() != '.') ) {
		file.pop_back();
	}

	return file;
}