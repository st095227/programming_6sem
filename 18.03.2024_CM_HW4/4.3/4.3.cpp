#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

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

double Left_rectangles(double A, double B, int m, double f(double))
{
	double h = (B - A) / m;
	double sum = 0;
	for (int i = 0; i < m; ++i)
	{
		sum += f(A + i * h);
	}
	return h * sum;
}

double Right_rectangles(double A, double B, int m, double f(double))
{
	double h = (B - A) / m;
	double sum = 0;
	for (int i = 0; i < m; ++i)
	{
		sum += f(A + (i + 1) * h);
	}
	return h * sum;
}

double Middle_rectangles(double A, double B, int m, double f(double))
{
	double h = (B - A) / m;
	double sum = 0;
	for (int i = 0; i < m; ++i)
	{
		sum += f(A + (i + 0.5) * h);
	}
	return h * sum;
}
double trapezoids(double A, double B, int m, double f(double))
{
	double h = (B - A) / m;
	double sum = 0;
	sum += f(A);
	for (int i = 1; i < m; ++i)
	{
		sum += 2 * f(A + i * h);
	}
	sum += f(B);
	return 0.5 * h * sum;
}

double SimpsonM(double A, double B, int m, double f(double))
{
	double h = (B - A) / m;
	double sum = 0;
	sum += f(A);
	for (int i = 1; i < m; i = i + 2)
	{
		sum += 4 * f(A + i * h);
	}
	for (int i = 2; i < m - 1; i = i + 2)
	{
		sum += 2 * f(A + i * h);
	}
	sum += f(B);
	return  (h / 3) * sum;
}

double Runge(double a, double b, double m, double d, double met(double, double, int, double(*)(double)), double f(double))
{
	return (pow(2, d + 1) * met(a, b, 2 * m, f) - met(a, b, m, f)) / (pow(2, d + 1) - 1);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Приближённое вычисление интеграла по составным квадратурным формулам." << endl;
	cout << "Введите границы промежутка интегрирования и количество промежутков разбиения: ";
	
	double a, b;
	int m, l;
	cin >> a >> b >> m;
	cout << endl;
	cout << "Введите параметр увеличения количества промежутков: ";
	cin >> l;
	cout << endl;
	//cout << "Введите способ:" << endl << "1.КФ левых прямоугольников" << endl << "2.КФ правых прямоугольников" << endl << "3.КФ средних прямоугольников" << endl;
	//cout << "4.КФ трапеций" << endl << "5.КФ Симпсона" << endl;
	//int n;
	//cin >> n;
	double v1, v2, v3, v4, v5;

		v1 = Left_rectangles(a, b, m*l, func);
		cout << "Вычисленное значение КФ левых прямоугольников J(h/l): " << v1 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v1 - F(a, b)) << endl << endl;
	
		v2 = Right_rectangles(a, b, m*l, func);
		cout << "Вычисленное значение КФ правых прямоугольников J(h/l): " << v2 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v2 - F(a, b)) << endl << endl;

		v3 = Middle_rectangles(a, b, m*l, func);
		cout << "Вычисленное значение КФ средних прямоугольников J(h/l): " << v3 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v3 - F(a, b)) << endl << endl;

		v4 = trapezoids(a, b, m*l, func);
		cout << "Вычисленное значение КФ трапеций J(h/l): " << v4 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v4 - F(a, b)) << endl << endl;
		
		v5 = SimpsonM(a, b, m*l, func);
		cout << "Вычисленное значение КФ Симпсона J(h/l): " << v5 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v5 - F(a, b)) << endl << endl;
	///////////////////////////////////////////////
		v1 = Runge(a, b, m, 0, Left_rectangles,func);
		cout << "Уточненное значение КФ левых прямоугольников J(h): " << v1 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v1 - F(a, b)) << endl << endl;

		v2 = Runge(a, b, m, 0, Right_rectangles, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h): " << v2 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v2 - F(a, b)) << endl << endl;

		v3 = Runge(a, b, m, 0, Middle_rectangles, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h): " << v3 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v3- F(a, b)) << endl << endl;

		v4 = Runge(a, b, m, 1, trapezoids, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h): " << v4 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v4 - F(a, b)) << endl << endl;

		v5 = Runge(a, b, m, 3, SimpsonM, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h): " << v5 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v5 - F(a, b)) << endl << endl;
		//////////
		v1 = Runge(a, b, m * l, 0, Left_rectangles, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h/l): " << v1 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v1 - F(a, b)) << endl << endl;

		v2 = Runge(a, b, m * l, 0, Right_rectangles, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h/l): " << v2 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v2 - F(a, b)) << endl << endl;

		v3 = Runge(a, b, m * l, 0, Middle_rectangles, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h/l): " << v3 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v3 - F(a, b)) << endl << endl;

		v4 = Runge(a, b, m * l, 1, trapezoids, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h/l): " << v4 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v4 - F(a, b)) << endl << endl;

		v5 = Runge(a, b, m * l, 3, SimpsonM, func);
		cout << "Уточненное значение КФ левых прямоугольников J(h/l): " << v5 << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v5 - F(a, b)) << endl << endl;
	return 0;
}