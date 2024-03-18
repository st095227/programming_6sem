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
double Dfunc1(double x)
{
	return (double)-1 / (x * x);
}
double Dfunc2(double x)
{
	return (double)2 / (x * x * x);
}
double Dfunc3(double x)
{
	return (double)-6 / (x * x * x * x);
}
double Dfunc4(double x)
{
	return (double)24 / (x * x * x * x * x);
}
double F(double a, double b)
{
	return log(abs(b / a));
}

static bool abs_compare(double a, double b)
{
	return (abs(a) < abs(b));
}

double M(double a, double b, double f(double))
{
	vector<double> arr;
	double h = (b - a) / 1000;
	for(int i =0; i<=1000; ++i) arr.push_back(f(a + i * h));
	vector<double>::iterator result;
	result = max_element(arr.begin(), arr.end(), abs_compare);
	return arr[distance(arr.begin(), result)];
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
	for (int i = 2; i < m-1; i = i + 2)
	{
		sum += 2 * f(A +  i * h);
	}
	sum += f(B);
	return  (h/3) * sum;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Приближённое вычисление интеграла по составным квадратурным формулам." << endl;
	cout << "Введите границы промежутка интегрирования и количество промежутков разбиения: ";
	double a, b;
	int m;
	cin >> a >> b >> m;
	cout << endl;
	cout << "Введите способ:" << endl << "1.КФ левых прямоугольников" << endl << "2.КФ правых прямоугольников" << endl << "3.КФ средних прямоугольников" << endl;
	cout << "4.КФ трапеций" << endl << "5.КФ Симпсона" << endl;
	int n;
	cin >> n;
	double v = 0;
	double x;
	switch (n)
	{
	case 1:
		x = 0.5 * (b - a) * ((b - a) / m) * M(a, b, Dfunc1);
		v = Left_rectangles(a, b, m, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		cout << "теоретическая погрешность: " << x << endl;
		break;
	case 2:
		x = 0.5 * (b - a) * ((b - a) / m) * M(a, b, Dfunc1);
		v = Right_rectangles(a, b, m, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		cout << "теоретическая погрешность: " << x << endl;
		break;
	case 3:
		x = ((double)1/24) * (b - a) * ((b - a) / m) * M(a, b, Dfunc1);
		v = Middle_rectangles(a, b, m, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		cout << "теоретическая погрешность: " << x << endl;
		break;
	case 4:
		x = ((double)1 / 12) * (b - a) * ((b - a) / m)* ((b - a) / m) * M(a, b, Dfunc2);
		v = trapezoids(a, b, m, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		cout << "теоретическая погрешность: " << x << endl;
		break;
	case 5:
		x = ((double)1 / 2880) * (b - a) * ((b - a) / m)* ((b - a) / m)* ((b - a) / m)* ((b - a) / m) * M(a, b, Dfunc4);
		v = SimpsonM(a, b, m, func);
		cout << "Вычисленное значение: " << v << endl;
		cout << "Абсолютная фактическая погрешность: " << abs(v - F(a, b)) << endl;
		cout << "теоретическая погрешность: " << x << endl;
		break;
	}
	return 0;
}
