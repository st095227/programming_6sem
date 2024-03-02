#include<iostream>
#include<cmath>
#include<vector>
#include <algorithm>

using namespace std;

double f(double x)
{
	return sin(x) - (x*x)/2;
}

vector<vector<double>> table(double A, double B, int N)
{
	vector<vector<double>> tab;
	for (int i = 0; i <= N; ++i)
	{
		vector<double> arr = { ((B - A) * i) / N,f(((B - A) * i) / N) };
		tab.push_back(arr);
	}
	return tab;
}

bool compare_abs(double a, double b, double x)
{
	return(abs(a - x) < abs(b - x));
}

vector<double> sort_near(vector<double> arr, double x)
{
	double temp;
	for (int i = 0; i < arr.size(); i++) {
		for (int j = 0; j < arr.size() - i - 1; j++) {
			if (compare_abs(arr[j + 1], arr[j], x))
			{
				temp = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = temp;
			}
		}
	}
	return arr;
}

vector<double> matr_to_vec(vector<vector<double>> matr)
{
	vector<double> arr;
	for (int i = 0; i < matr.size(); ++i)
	{
		arr.push_back(matr[i][0]);
	}
	return arr;
}

double Lagrange(double x, vector<vector<double>> tab, int n)
{
	double sum = 0;
	vector<double> arr_sort = matr_to_vec(tab);
	arr_sort = sort_near(arr_sort, x);
	for (int k = 0; k <= n; ++k)
	{
		double temp = 1;
		double temp1 = 1;
		for (int i = 0; i < k; ++i)
		{
			temp = temp * (x - arr_sort[i]);
			temp1 = temp1 * (arr_sort[k] - arr_sort[i]);
		}
		for (int i = k+1; i <= n; ++i)
		{
			temp = temp * (x - arr_sort[i]);
			temp1 = temp1 * (arr_sort[k] - arr_sort[i]);
		}
		sum += temp * f(arr_sort[k]) / temp1;
		
	}
	return sum;
}

double Newton(double x, vector<vector<double>> tab,int n)
{
	vector<double> arr_sort = matr_to_vec(tab);
	arr_sort = sort_near(arr_sort, x);
	double sum = f(arr_sort[0]);
	for (int i = 1; i <= n; ++i) 
	{
		double F = 0;
		for (int j = 0; j <= i; ++j) 
		{   //следующее слагаемое полинома
			double den = 1;//считаем знаменатель разделенной разности
			for (int k = 0; k <= i; ++k)
				if (k != j)
					den *= (arr_sort[j] - arr_sort[k]);
			//считаем разделенную разность
			F += f(arr_sort[j]) / den;
		}

		//домножаем разделенную разность на скобки (x-x[0])...(x-x[i-1])
		for (int k = 0; k < i; ++k)
			F *= (x - arr_sort[k]);
		sum += F;//полином
	}
	return sum;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Задача алгебраического интерполирования. Вариант 1." << endl;
	cout << "Введите число значений в таблице: ";
	int m;
	cin >> m;
	cout << endl << "Введите концы отрезка, из которого выбираются узлы: ";
	double a, b;
	cin >> a >> b;
	cout << endl;
	vector<vector<double>> tab = table(a, b, m-1);

	cout << "X(i)" << ' ' << "f(X)" << endl;
	for (int i = 0; i < tab.size(); ++i)
	{
		cout << tab[i][0] << ";  " << tab[i][1] << endl;
	}

	double x;
	cout << endl << "Введите точку интерполирования: ";
	cin >> x;
	int n;
	while (1)
	{
		cout << endl << "Введите степень интерполяционного многочлена: ";
		cin >> n;

		if (n == -1) break;

		vector<double> arr = matr_to_vec(tab);

		cout << "набор узлов интерполирования:" << endl;
		vector<double> arrsort = sort_near(arr, x);
		for (int i = 0; i < n; ++i)
		{
			cout << arrsort[i] << " ";
		}
		cout << endl;
		cout << endl << "lagrange:" << Lagrange(x, tab, n) << endl;
		cout << "Значение абсолютной фактической погрешности для формы Лагранжа: " << abs(f(x) - Lagrange(x, tab, n)) << endl << endl;
		cout << endl << "Newton:" << Newton(x, tab, n) << endl;
		cout << "Значение абсолютной фактической погрешности для формы Ньютона: " << abs(f(x) - Newton(x, tab, n)) << endl << endl;
		cout << "для выхода введите -1." << endl << endl;
	}
	return 0;
}