#include<iostream>
#include<cmath>
#include<vector>
using namespace std;

double f(double x)
{
	return x - 10 * sin(x);
}

vector<double> separate_roots(double A, double B, int N)
{
	double H = (B - A) / N;
	int counter = 0;
	double x1 = A;
	double x2 = x1 + H;
	double y1 = f(x1);
	vector<double> arr;
	while (x2 <= B)
	{
		double y2 = f(x2);
		if (y1 * y2 <= 0)
		{
			++counter;
			arr.push_back(x1);
			arr.push_back(x2);
		}
		x1 = x2;
		x2 = x1 + H;
		y1 = y2;
	}
	arr.push_back(counter);
	return arr;
}

int bisection(double A, double B, double eps)
{
	double a = A;
	double b = B;
	int n = 0;
	while (b - a > 2 * eps)
	{
		double c = (a + b) / 2;
		if (f(a) * f(c) <= 0)
		{
			b = c;
		}
		else
		{
			a = c;
		}
		++n;
	} 
	cout << "Начальное приближение: " << (A + B) / 2 << endl;
	cout << " Корень на промежутке [" << A << ";" << B << "]: " << (b + a) / 2 << endl;
	cout << " |Xm-X(m-1)|: " << (b - a) / 2 << endl;
	cout << " Величина невязки |f(x)-0|: " << abs(f((b + a) / 2) - 0) << endl << endl;
	return n;
}

double deriv(double x, double f1(double x))
{
	const double h = 1e-6;
	return (f1(x + h) - f1(x - h)) / (2.0 * h);
}
double deriv2(double x, double f1(double x))
{
	const double h = 1e-6;
	return (f1(x + h) - 2.0 * f1(x) + f1(x - h)) / (h * h);
}

int newtone(double A, double B, double eps)
{
	double a = A;
	double b = B;
	int n = 0;
	double x0 = (b - a) / 2;
	if (f(A) * deriv2(A, f) > 0)
	{
		x0 = A;
	}
	else if (f(B) * deriv2(B, f) > 0)
	{
		x0 = B;
	}
	else
	{
		cout << "не удалось найти решение на промежутке" << "[" << A << "; " << B << "]" << endl << endl;
		return 0;
	}
	cout << "Начальное приближение: " << x0 << endl;
	++n;
	double x = x0 - (f(x0) / deriv(x0, f));
	while (abs(x - x0) > eps)
	{
		++n;
		x0 = x;
		x = x0 - (f(x0) / deriv(x0, f));
	}
	cout << " Корень на промежутке [" << A << ";" << B << "]: " << x << endl;
	cout << " |Xm-X(m-1)|: " << abs(x - x0) << endl;
	cout << " Величина невязки |f(x)-0|: " << abs(f(x) - 0) << endl << endl;
	return n;
}

int mod_newtone(double A, double B, double eps)
{
	double a = A;
	double b = B;
	int n = 0;
	double x0 = (b-a)/2;
	if (f(A) * deriv2(A, f) > 0)
	{
		x0 = A;
	}
	else if (f(B) * deriv2(B, f) > 0)
	{
		x0 = B;
	}
	else
	{
		cout << "не удалось найти решение на промежутке" << "[" << A << "; " << B << "]" << endl << endl;
		return 0;
	}
	cout << "Начальное приближение: " << x0 << endl;
	++n;
	double xk = x0;
	double x = x0 - (f(x0) / deriv(x0, f));
	while (abs(x - xk) > eps)
	{
		++n;
		xk = x;
		x = xk - (f(xk) / deriv(x0, f));
	}
	cout << " Корень на промежутке [" << A << ";" << B << "]: " << x << endl;
	cout << " |Xm-X(m-1)|: " << abs(x - xk) << endl;
	cout << " Величина невязки |f(x)-0|: " << abs(f(x) - 0) << endl << endl;
	return n;
}

int secant(double A, double B, double eps)
{
	double a = A;
	double b = B;
	int n = 0;
	double x0 = 0;
	double x1 = 0;
	x0 = a;
	x1 = b;
	cout << "Начальное приближение: " << a << " " << b << endl;
	++n;
	double xk = x1- (f(x1)/(f(x1)-f(x0)))*(x1-x0);
	while (abs(xk - x1) > eps)
	{
		++n;
		x0 = x1;
		x1 = xk;
		xk = x1 - (f(x1) / (f(x1) - f(x0))) * (x1 - x0);
	}
	cout << " Корень на промежутке [" << A << ";" << B << "]: " << xk << endl;
	cout << " |Xm-X(m-1)|: " << abs(xk - x1) << endl;
	cout << " Величина невязки |f(x)-0|: " << abs(f(xk) - 0) << endl << endl;
	return n;
}

int main()
{
	setlocale(LC_ALL, "Russian");

	double A = -5;
	double B = 3;
	double eps = 1e-6;
	int N = 100;

	cout << "Решение нелинейных уравнений с исходными параметрами." << endl;
	cout << "Исходные параметры:" << endl << "Уравнение x - 10 * sin(x) = 0." << endl;
	cout << "Корни на промежутке [" << A << ";" << B << "] c погрешностью " << eps << "." << endl;
	cout << "Шаг табулирования промежутка: " << (B - A) / N << endl << endl;


	vector<double> segments = separate_roots(A, B, N);
	cout << "Количество знакопеременных промежутков: " << segments[segments.size() - 1] << endl;
	for (int i = 0; i < segments.size()-1; ++i)
	{
		cout << '[' << segments[i] << " , " << segments[i + 1] << ']' << endl;
		++i;
	}

	cout << endl << "Определение корней методом БИССЕКЦИИ" << endl;
	int k = 0;
	for (int i = 0; i < segments.size()-1; i=i+2)
	{
		k += bisection(segments[i], segments[i+1], eps);
	}
	cout << " Количество шагов для определения всех корней методом бисекции: " << k << endl << endl;

	cout << endl << "Определение корней методом НЬЮТОНА" << endl;
	k = 0;
	for (int i = 0; i < segments.size() - 1; i = i + 2)
	{
		k += newtone(segments[i], segments[i + 1], eps);
	}
	cout << " Количество шагов для определения всех корней методом Ньютона: " << k << endl << endl;

	cout << endl << "Определение корней модифицированным методом НЬЮТОНА" << endl;
	k = 0;
	for (int i = 0; i < segments.size() - 1; i = i + 2)
	{
		k += mod_newtone(segments[i], segments[i + 1], eps);
	}
	cout << " Количество шагов для определения всех корней модифицированным методом Ньютона: " << k << endl << endl;

	cout << endl << "Определение корней методом секущих" << endl;
	k = 0;
	for (int i = 0; i < segments.size() - 1; i = i + 2)
	{
		k += secant(segments[i], segments[i + 1], eps);
	}
	cout << " Количество шагов для определения всех корней методом секущих: " << k << endl << endl;
	return 0;
}

