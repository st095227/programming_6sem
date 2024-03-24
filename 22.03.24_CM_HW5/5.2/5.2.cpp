#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const double Integ_const = 1.605412976803;
const double Integ_const2 = 2.403939430634;
double func(double x)
{
	return (sin(x) / x);
}
double func5(double x)
{
	return pow(x, 5);
}
double Legendre(double n, double x)
{
	if (n == 0)
	{
		return 1;
	}
	else if (n == 1)
	{
		return x;
	}
	else
	{
		return (((2 * n - 1) / n) * Legendre(n - 1, x)*x - (((n - 1) / n) * Legendre(n - 2, x)));
	}
}

double secantLeg(double a, double b, double eps, int n)
{

	double x0 = a;
	double x1 = b;
	double xk = x1 - (Legendre(n, x1) / (Legendre(n, x1) - Legendre(n, x0))) * (x1 - x0);
	while (abs(xk - x1) > eps)
	{
		x0 = x1;
		x1 = xk;
		xk = x1 - (Legendre(n, x1) / (Legendre(n, x1) - Legendre(n, x0))) * (x1 - x0);
	}
	return xk;
}

vector<double> separate_roots(int N, int n)
{
	double H = (double)2 / N;
	double x1 = -1;
	double x2 = x1 + H;
	double y1 = Legendre(n, x1);
	vector<double> arr;
	while (x2 <= 1)
	{
		double y2 = Legendre(n, x2);
		if (y1 * y2 <= 0)
		{
			arr.push_back(x1);
			arr.push_back(x2);
		}
		x1 = x2;
		x2 = x1 + H;
		y1 = y2;
	}
	return arr;
}

double KF_Gauss(double f(double), double a, double b, double n)
{
	vector<double> segm;
	vector<double> roots;
	vector<double> coeff;
	segm = separate_roots(1000, n);
	for (int i = 0; i < segm.size(); i = i + 2)
	{
		roots.push_back(secantLeg(segm[i], segm[i + 1], 1e-12, n));
	}
	for (int i = 0; i < roots.size(); ++i)
	{
		coeff.push_back(((2 * (1 - roots[i] * roots[i])) / (n * n * Legendre(n - 1, roots[i]) * Legendre(n - 1, roots[i]))));
	}
	cout << "x(i)   a(i)" << endl;
	for (int i = 0; i < roots.size(); ++i)
	{
		cout << roots[i] << ' ' << coeff[i] << endl;
	}
	double sum = 0;
	for (int i = 0; i < coeff.size(); ++i)
	{
		sum += coeff[i] * f(0.5*(a + b + (b - a) * roots[i]));
	}
	sum *= (0.5 * (b - a));
	return sum;
}
double MELER(double f(double), double N)
{
	cout << "Узлы(корни многочлена Чебышева): " << endl;
	vector<double> roots;
	for (int i = 1; i <= N; ++i)
	{
		roots.push_back(cos(3.14 * (2 * i - 1) / (2 * N)));
	}
	for (int i = 0; i < roots.size(); ++i)
	{
		cout << roots[i] << endl;
	}
	double sum = 0;
	for (int i = 0; i < roots.size(); ++i)
	{
		sum += f(roots[i]);
	}
	return sum *= (3.14 / N);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "КФ Гаусса." << endl;
	cout << "Введите границы промежутка интегрирования: ";
	double a, b;
	cin >> a >> b;
	cout << endl;
	cout << "Точное значение: " << Integ_const << endl << endl;
	double ans = KF_Gauss(func, a, b, 5);
	cout << endl << "Найденное значение: " << ans;
	cout << endl << "Модуль невязки: " << abs(ans - Integ_const);
	cout << endl << endl;
	cout << "Вычисление интегралов при помощи КФ Мелера." << endl;
	double N1, N2, N3;
	cout << "Введите количество узлов (3 значения): ";
	cin >> N1 >> N2 >> N3;
	cout << endl;
	cout << "значение интеграла при N1 узлах: " << MELER(cos, N1) << endl << endl;
	cout << "значение интеграла при N2 узлах: " << MELER(cos, N2) << endl << endl;
	cout << "значение интеграла при N3 узлах: " << MELER(cos, N3) << endl << endl;
	return 0;
}