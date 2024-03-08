#define _USE_MATH_DEFINES
#include<iostream>
#include<cmath>
#include<vector>
#include <algorithm>

using namespace std;

double f(double x)
{
	return pow(M_E,3*x);
}
double der1(double x)
{
	return 3 * pow(M_E, 3 * x);
}
double der2(double x)
{
	return 9 * pow(M_E, 3 * x);
}
vector<vector<double>> table(double A, double h, int N)
{
	vector<vector<double>> tab;
	for (int i = 0; i <= N; ++i)
	{
		vector<double> arr = { A+i*h , f(A + i * h) };
		tab.push_back(arr);
	}
	return tab;
}

vector<double> derND(double a, double h, double m)
{
	vector<double> ans;
	double fd = 0;
	fd = (-3 * f(a) + 4 * f(a + h) - f(a + 2 * h)) / (2 * h);
	ans.push_back(fd);
	for (int i = 2; i < m; ++i)
	{
		ans.push_back((f(a + i * h) - f(a + (i - 2) * h)) / (2 * h));
	}
	ans.push_back((3 * f(a + h * (m - 1)) - 4 * f(a + h * (m - 2)) + f(a + h * (m - 3))) / (2 * h));
	return ans;
}

vector<double> err1(vector<double> dev, double a, double h)
{
	vector<double> ans;
	for (int i = 0; i < dev.size(); ++i) ans.push_back(abs(dev[i]-der1(a+i*h)));
	return ans;
}
vector<double> derND2(vector<double> devD, double h)
{
	vector<double> ans;
	double fd = 0;
	fd = (-3 * devD[0] + 4 * devD[1] - devD[2]) / (2 * h);
	ans.push_back(fd);
	for (int i = 1; i < devD.size()-1; ++i)
	{
		ans.push_back((devD[i+1] - devD[i-1]) / (2 * h));
	}
	ans.push_back((3 * devD[devD.size() - 1] - 4 * devD[devD.size() - 2] + devD[devD.size() - 3]) / (2 * h));
	return ans;
}

vector<double> err2(vector<double> dev, double a, double h)
{
	vector<double> ans;
	for (int i = 0; i < dev.size(); ++i) 
{
		ans.push_back(abs(dev[i] - der2(a + i * h)));
	}
	return ans;
}


int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Задача нахождения производных по формулам численного интегрирования." << endl;
	cout << "Введите число значений в таблице: ";
	int m;
	cin >> m;
	cout << endl << "Введите начальную точку:";
	double a;
	cin >> a;
	cout << endl << "Введите шаг:";
	double h;
	cin >> h;
	vector<vector<double>> tab = table(a, h, m - 1);
	cout << "X(i)" << ' ' << "f(x)" << endl;
	for (int i = 0; i < tab.size(); ++i)
	{
		cout << tab[i][0] << ";  " << tab[i][1] << endl;
	}
	cout << endl;
	vector<double> deriv1 = derND(a, h, m);
	vector<double> deriv2 = derND2(deriv1, h);
	vector<double> error1 = err1(deriv1, a, h);
	vector<double> error2 = err2(deriv2, a, h);
	for (int i = 0; i < m; ++i)
	{
		cout << deriv1[i] << ";  " << error1[i] << "; " << deriv2[i] << ";  " << error2[i] <<endl;
	}
	return 0;
}