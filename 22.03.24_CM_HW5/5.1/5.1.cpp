#include <iostream>
#include <cmath>
#include <vector>

using namespace std;
const double Integ_const = 0.364221932032; // sin(x)*sqrt(x)

double func(double x)
{
	return sin(x);
}
double f1(double x)
{
	return x;
}
double f2(double x)
{
	return 1;
}
double poly(double a, double b, double c, double d, double x)
{
	return x * x * x * d + x * x * c + x * b + a;
}

double func1(double x)
{
	return poly(5, 0, 0, 0, x);
}
double func2(double x)
{
	return poly(5, 2, 0, 0, x);
}

double IKF_2node(double f(double), double c, double d, double x1, double x2)
{
	cout << "mom0: " << ((double)2 / 3) * (pow(d, 1.5) - pow(c, 1.5)) << endl;
	cout << "mom1: " << ((double)2 / 5) * (pow(d, 2.5) - pow(c, 2.5)) << endl;
	cout << "mom2: " << ((double)2 / 7) * (pow(d, 3.5) - pow(c, 3.5)) << endl;
	double A2 = (6 * pow(d, 2.5) - 10 * x1 * pow(d, 1.5) - 6 * pow(c, 2.5) + 10 * x1 * pow(c, 1.5)) / (15 * x2 - 15 * x1);
	double A1 = (6 * pow(d, 2.5) - 10 * x2 * pow(d, 1.5) - 6 * pow(c, 2.5) + 10 * x2 * pow(c, 1.5)) / (15 * x1 - 15 * x2);
	cout << "A1: " << A1 << endl;
	cout << "A2: " << A2 << endl;
	cout << "Приближенное значение ИКФ по 2 узлам: " << A1 * f(x1) + A2 * f(x2);
	return A1 * f(x1) + A2 * f(x2);
}
 
double KF_GT_2node(double f(double), double c, double d)
{
	double m0 = ((double)2 / 3) * (pow(d, 1.5) - pow(c, 1.5));
	double m1 = ((double)2 / 5) * (pow(d, 2.5) - pow(c, 2.5));
	double m2 = ((double)2 / 7) * (pow(d, 3.5) - pow(c, 3.5));
	double m3 = ((double)2 / 9) * (pow(d, 4.5) - pow(c, 4.5));
	cout << "mom0: " << m0 << endl;
	cout << "mom1: " << m1 << endl;
	cout << "mom2: " << m2 << endl;
	cout << "mom3: " << m3 << endl;

	double a1 = ((m0 * m3) - (m2 * m1)) / ((m1 * m1) - m2 * m0);
	double a2 = ((m2 * m2) - (m3 * m1)) / ((m1 * m1) - m2 * m0);
	
	double x1 = 0.5 * (-a1 - sqrt((a1 * a1) - 4 * a2));
	double x2 = 0.5 * (-a1 + sqrt((a1 * a1) - 4 * a2));
	cout << "x1: " << x1 << endl;
	cout << "x2: " << x2 << endl;

	double A1 = ((double)1 / (x1 - x2)) * (m1 - x2 * m0);
	double A2 = ((double)1 / (x2 - x1)) * (m1 - x1 * m0);
	cout << "Приближенное значение ИКФ Типа Гаусса: " << A1 * f(x1) + A2 * f(x2);
	return A1 * f(x1) + A2 * f(x2);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Приближённое вычисление интегралов при помощи КФ НАСТ." << endl;
	cout << "Введите границы промежутка интегрирования: ";
	double a, b;
	cin >> a >> b;
	cout << endl;
	cout << "Точное значение: " << Integ_const << endl << endl;
	cout << "ИКФ по 2 узлам" << endl;
	IKF_2node(func, a, b, a, b);
	cout << endl;
	KF_GT_2node(func, a, b);
	return 0;
}