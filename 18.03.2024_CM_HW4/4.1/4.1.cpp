#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double poly(double a, double b, double c, double d, double x)
{
	return x * x * x * d + x * x * c + x * b + a;
}

double func1(double x)
{
	return poly(5,0,0,0,x);
}
double func2(double x)
{
	return poly(5, 2, 0, 0, x);
}
double func3(double x)
{
	return poly(5, 2, 3, 0, x);
}
double func4(double x)
{
	return poly(5, 2, 3, 4, x);
}

double F1(double a, double b)
{
	return 5 * (b - a);
}
double F2(double a, double b)
{
	return 5 * (b - a) + (b * b - a * a);
}
double F3(double a, double b)
{
	return 5 * (b - a) + (b * b - a * a) + (b * b * b - a * a * a);
}
double F4(double a, double b)
{
	return 5 * (b - a) + (b * b - a * a) + (b * b * b - a * a * a) + (b * b * b * b - a * a * a * a);
}

double func(double x)
{
	return (double)1 / x;
}
double F(double a, double b)
{
	return log(abs(b / a));
}

double Left_rectangle(double a, double b, double f(double))
{
	return (b - a) * f(a);
}
double Right_rectangle(double a, double b, double f(double))
{
	return (b - a) * f(b);
}
double Middle_rectangle(double a, double b, double f(double))
{
	return (b - a) * f((b + a) / 2);
}
double trapezoid(double a, double b, double f(double))
{
	return ((b - a) / 2) * (f(a) + f(b));
}
double Simpson(double a, double b, double f(double))
{
	return ((b - a) / 6) * (f(a) + 4 * f((b + a) / 2) + f(b));
}
double KF38 (double a, double b, double f(double))
{
	double h = (b - a) / 3;
	return (b - a) * ((0.125) * f(a) + (0.375) * f(a + h) + (0.375) * f(a + 2 * h) + (0.125) * f(b));
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Приближённое вычисление интеграла по квадратурным формулам." << endl;
	cout << "Введите границы промежутка интегрирования: ";
	double a, b;
	cin >> a >> b;
	cout << endl;
	cout << "Введите способ:" << endl << "1.КФ левого прямоугольника" << endl << "2.КФ правого прямоугольника" << endl << "3.КФ среднего прямоугольника" << endl;
	cout << "4.КФ трапеции" << endl << "5.КФ Симпсона" << endl << "6.КФ 3/8" << endl;
	int n;
	cin >> n;
	double v = 0;
	switch (n)
	{
	case 1:
		v = Left_rectangle(a, b, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		break;
	case 2:
		v = Right_rectangle(a, b, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		break;
	case 3:
		v = Middle_rectangle(a, b, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		break;
	case 4:
		v = trapezoid(a, b, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		break;
	case 5:
		v = Simpson(a, b, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		break;
	case 6:
		v = KF38(a, b, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		break;
	}
	return 0;
}