#include "ResFunc.h"

using namespace std;

ResFunc::ResFunc()
{
}

int col_count, row_count;
double A, B, C, D, E, F;

double poly(double p, double T)
{
	return A+B*p+C*T+D*p*p+E*T*T+F*p*T;
}

Matrix Read(string s, int m, int n)
{
	//create array a with m columns and n rows
	Matrix a(m, vector<double>(n));

	col_count = m;
	row_count = n;

	// create file-stream for @table file - file with f(p,T) values in tabular form
	ifstream in(s);
	if (in.is_open())
	{
		//counting amount of numbers in file
		int count = 0;
		int temp;

		while(!in.eof())
		{
			in >> temp; //read table
			count++; //amount of numbers
		}

		//TAB-counter
		//go to beginning of file
		in.seekg(0, ios::beg);
		in.clear();

		int tab_count = 0; // TAB-count = 0
		char symbol;
		while(!in.eof())
		{
			in.get(symbol); // read current symbol
			if (symbol == '	') tab_count++;
			if (symbol == '\n') break;
		}

		//go to beginning of file again
		in.seekg(0, ios::beg);
		in.clear();

		col_count = count / (tab_count + 1); // col number
		row_count = tab_count + 1; // row number

		//read table from file
		for (int i = 0; i<row_count; i++)
			for (int j = 0; j<col_count; j++)
				in >> a[i][j];
		
		/* destructor
		for (int i = 0; i<n; i++) delete [] a[i];
		delete[] a;
		*/
		in.close();
	}
	else 
	{
		cout << "Problems with open table file";
	}

	//...///
}

Matrix b(col_count, vector<double>(row_count));

double value(const Vector& X) 
{
	double pressures[10] = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200};
	double temperatures[6] = {293.15, 303.15, 313.15, 321.65, 333.15, 343.15};

	int ii, jj, c, sc;

	A = X[0];
	B = X[1];
	C = X[2];
	D = X[3];
	E = X[4];
	F = X[5];

	for(double tt:temperatures)
	{
		for(double pp:pressures) 
		{
			b[ii][jj] = poly(pp, tt);
			jj++;
		}
		ii++;
	}

	for (int i = 0; i<row_count; i++)
		for (int j = 0; j<col_count; j++)
		{
			c = a[i][j]-b[i][j];
			sc += c*c;
		}
		return sc;
}

ResFunc::~ResFunc(void)
{
}