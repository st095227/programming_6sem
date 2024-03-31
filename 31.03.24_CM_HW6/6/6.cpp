#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double y(double x)
{
	return (pow(M_E, 2 * x) - 1) / (pow(M_E, 2 * x) + 1);
}

double f(double X, double Y)
{
	return -pow(Y, 2) + 1;
}

void out_accurate_solution_table(double step, int qt, double x0)
{
	for (int i = -2; i <= qt; ++i)
	{
		cout << x0 + i * step << ' ' << y(x0 + i * step) << endl;
	}
}

vector<vector<double>> Euler1(double step, int qt, double x0, double y0, double func(double, double))
{
	vector<vector<double>> ans;
	ans.push_back({x0,y0});
	double temp_y_12 = 0;
	double temp_yk1 = y0;
	for (int i = 0; i < qt; ++i)
	{
		temp_y_12 = temp_yk1 + (step / 2) * func(x0 + i * step, temp_yk1);
		temp_yk1 = temp_yk1 + step * func(x0 + i * step + (step / 2), temp_y_12);
		ans.push_back({ x0 + (i+1) * (step), temp_yk1 });
	}
	return ans;
}

vector<vector<double>> Euler2(double step, int qt, double x0, double y0, double func(double, double))
{
	vector<vector<double>> ans;
	ans.push_back({ x0,y0 });
	double temp_Yk1 = 0;
	double temp_yk1 = y0;
	for (int i = 0; i < qt; ++i)
	{
		temp_Yk1 = temp_yk1 + (step) * func(x0 + i * step, temp_yk1);
		temp_yk1 = temp_yk1 + (step / 2) * (func(x0 + i * step, temp_yk1) + f(x0 + (i + 1) * step, temp_Yk1));
		ans.push_back({ x0 + (i + 1) * (step), temp_yk1 });
	}
	return ans;
}

vector<vector<double>> Runge_Kutta_method(double step, int qt, double x0, double y0, double func(double, double))
{
	vector<vector<double>> ans;
	ans.push_back({ x0,y0 });
	double k1, k2, k3, k4;
	double yk = y0;
	for (int i = 0; i < qt; ++i)
	{
		k1 = step * func(x0 + i * step, yk);
		k2 = step * func(x0 + i * step + (step / 2), yk + (k1 / 2));
		k3 = step * func(x0 + i * step + (step / 2), yk + (k2 / 2));
		k4 = step * func(x0 + i * step + step, yk + k3);
		yk = yk + ((double)1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
		ans.push_back({ x0 + (i + 1) * (step), yk });
	}
	return ans;
}

double der1(double X, double Y)
{
	return f(X, Y);
}
double der2(double X, double Y)
{
	return -2 * Y * der1(X, Y);
}
double der3(double X, double Y)
{
	return -2 * Y * der2(X, Y) - 2 * der1(X, Y) * der1(X, Y);
}
double der4(double X, double Y)
{
	return -2 * (Y * der3(X, Y) + der1(X, Y) * der2(X, Y) + 2 * der1(X, Y) * der2(X, Y));
}
double der5(double X, double Y)
{
	return -2 * (Y * der4(X, Y) + 4 * der1(X, Y) * der3(X, Y) + 3 * der2(X, Y) * der2(X, Y));
}

vector<vector<double>> Taylor(double step, int qt, double x0, double y0, double func(double, double))
{
	vector<vector<double>> ans;
	double yk;
	for (int i = -2; i <= qt; ++i)
	{
		yk = y0 + (der1(x0, y0) * (x0 + i * step - x0)) + ((der2(x0, y0) / 2) * pow((x0 + i * step - x0), 2)) + ((der3(x0, y0) / 6) * pow((x0 + i * step - x0), 3)) + ((der4(x0, y0) / 24) * pow((x0 + i * step - x0), 4)) + ((der5(x0, y0) / 120) * pow((x0 + i * step - x0), 5));
		ans.push_back({ x0 + i * step, yk });
	}
	return ans;
}
double Finite_difference(vector<double> yk, int n, int k)
{
	if (n == 1)
	{
		return yk[k + 1] - yk[k];
	}
	else
	{
		return Finite_difference(yk, n - 1, k + 1) - Finite_difference(yk, n - 1, k);
	}
}

vector<vector<double>> Adams(double step, int qt, double x0, double y0, double func(double, double))
{
	vector<vector<double>> rez;
	rez = Taylor(step, 5, x0, y0, func);
	vector<vector<double>> ans;
	vector<double> yk;
	for (int i = 0; i < rez.size(); ++i)
	{
		yk.push_back(rez[i][1]);
	}
	double yk1 = yk[4];
	vector<double> qk;
	for (int i = 0; i <= 4; ++i)
	{
		qk.push_back(step * func(x0 + i * step, yk[i]));
	}

	for (int i = 4; i < qt+2; ++i)
	{
		yk1 = yk1 + qk[i] + ((double)1 / 2) * Finite_difference(qk, 1, i - 1) + ((double)5 / 12) * Finite_difference(qk, 2, i - 2) + ((double)3 / 8) * Finite_difference(qk, 3, i - 3) + ((double)251 / 720) * Finite_difference(qk, 4, i - 4);
		qk.push_back(step * func(x0 + (i+1) * step, yk1));
		ans.push_back({ x0 + (i-1) * step, yk1 });
	}
	return ans;
}



int main()
{

	double h = 0.1;
	double N = 10;
	cout << "Accurate solution:" << endl;
	out_accurate_solution_table(h, N, 0);

	vector<vector<double>> rez;
	cout << endl << "Taylor:" << endl;
	rez = Taylor(h, N, 0, 0, f);
	for (int i = 0; i < rez.size(); ++i)
	{
		cout << rez[i][0] << ' ' << rez[i][1] << endl;
	}
	cout << endl;
	cout << "Absolute errors of the Taylor method:" << endl;
	for (int i = 0; i < rez.size(); ++i)
	{
		cout << abs(rez[i][1]-y(rez[i][0])) << endl;
	}
	cout << endl;

	rez = Adams(h, N, 0, 0, f);
	cout << "Adams method:" << endl;
	for (int i = 0; i < rez.size(); ++i)
	{
		cout << rez[i][0] << ' ' << rez[i][1] << endl;
	}
	cout << "last value error:" << abs(rez[rez.size() - 1][1] - y(rez[rez.size() - 1][0])) << endl << endl;
	
	rez = Runge_Kutta_method(h, N, 0, 0, f);
	cout << "Runge Kutta method:" << endl;
	for (int i = 0; i < rez.size(); ++i)
	{
		cout << rez[i][0] << ' ' << rez[i][1] << endl;
	}
	cout << "last value error:" << abs(rez[rez.size() - 1][1] - y(rez[rez.size() - 1][0])) << endl << endl;


	rez = Euler1(0.1, 10, 0, 0, f);
	cout << "Euler I method:" << endl;
	for (int i = 0; i < rez.size(); ++i)
	{
		cout << rez[i][0] << ' ' << rez[i][1] << endl;
	}
	cout << "last value error:" << abs(rez[rez.size() - 1][1] - y(rez[rez.size() - 1][0])) << endl << endl;


	rez = Euler2(0.1, 10, 0, 0, f);
	cout << "Euler II method:" << endl;
	for (int i = 0; i < rez.size(); ++i)
	{
		cout << rez[i][0] << ' ' << rez[i][1] << endl;
	}
	cout << "last value error:" << abs(rez[rez.size() - 1][1] - y(rez[rez.size() - 1][0])) << endl;

	return 0;
}