#include<iostream>
#include<cmath>
#include<vector>
#include <algorithm>

using namespace std;

double f(double x)
{
	return sin(x) - (x * x) / 2;
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

bool compare_abs(vector<double> a, vector<double> b, double x)
{
	return(abs(a[0] - x) < abs(b[0] - x));
}

/*vector<double> sort_near(vector<double> arr, double x)
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
}*/

vector<vector<double>> sort_near(vector<vector<double>> arr, double x)
{
	vector<double> temp;
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

/*double Lagrange(double x, vector<vector<double>> tab, int n)
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
		for (int i = k + 1; i <= n; ++i)
		{
			temp = temp * (x - arr_sort[i]);
			temp1 = temp1 * (arr_sort[k] - arr_sort[i]);
		}
		sum += temp * f(arr_sort[k]) / temp1;

	}
	return sum;
}*/

double Newton(double x, vector<vector<double>> tab, int n)
{
	vector<vector<double>> arr_sort = sort_near(tab, x);
	//for (int i = 0; i < arr_sort.size(); ++i)
	//{
	//	cout << arr_sort[i][0] << " " << arr_sort[i][1] << endl;
	//}
	double sum = arr_sort[0][1];
	for (int i = 1; i <= n; ++i)
	{
		double F = 0;
		for (int j = 0; j <= i; ++j)
		{   //следующее слагаемое полинома
			double den = 1;//считаем знаменатель разделенной разности
			for (int k = 0; k <= i; ++k)
				if (k != j)
					den *= (arr_sort[j][0] - arr_sort[k][0]);
			//считаем разделенную разность
			F += arr_sort[j][1] / den;
		}

		//домножаем разделенную разность на скобки (x-x[0])...(x-x[i-1])
		for (int k = 0; k < i; ++k)
			F *= (x - arr_sort[k][0]);
		sum += F;//полином
	}
	return sum;
}

vector<vector<double>> reverse_table(double A, double B, int N)
{
	vector<vector<double>> tab;
	for (int i = 0; i <= N; ++i)
	{
		vector<double> arr = { f(((B - A) * i) / N),((B - A) * i) / N };
		tab.push_back(arr);
	}
	return tab;
}

double backward_interpolation_w1(double x, vector<vector<double>>  rtab, int n)
{
	return Newton(x, rtab, n);
}

vector<vector<double>> tabPol(double A, double B, int N, vector<vector<double>> table1, int k, double F)
{
	vector<vector<double>> tab;
	for (int i = 0; i <= N; ++i)
	{
		vector<double> arr = { ((B - A) * i) / N, Newton(((B - A) * i) / N, table1, k)-F};
		tab.push_back(arr);
	}
	return tab;
}

vector<double> separate_roots_for_Pol(vector<vector<double>> Pol_tab)
{
	double H = Pol_tab.size();
	vector<double> arr;
	for(int i=0; i< Pol_tab.size()-1; ++i)
	{
		
		if ((Pol_tab[i][1]) * (Pol_tab[i+1][1]) <= 0)
		{
			arr.push_back(Pol_tab[i][0]);
			arr.push_back(Pol_tab[i+1][0]);
		}
		
	}
	return arr;
}

double bisection(double A, double B, double eps, vector<vector<double>> table1, double F)
{
	double a = A;
	double b = B;

	while (b - a > 2 * eps)
	{
		double c = (a + b) / 2;
		if ((Newton(a, table1, 10) - F) * (Newton(c, table1, 10) - F) <= 0)
		{
			b = c;
		}
		else
		{
			a = c;
		}
	}
	//cout << "Ќачальное приближение: " << (A + B) / 2 << endl;
	//cout << "  орень на промежутке [" << A << ";" << B << "]: " << (b + a) / 2 << endl;
	//cout << " |Xm-X(m-1)|: " << (b - a) / 2 << endl;
	//cout << " ¬еличина нев€зки |f(x)-0|: " << abs(f((b + a) / 2) - 0) << endl << endl;
	return (b + a) / 2;
}

vector<double> backward_interpolation_w2(double F, vector<vector<double>> tab, int n, double eps)
{
	double x = 0;
	vector<double> intervals = separate_roots_for_Pol(tabPol(tab[0][0], tab[tab.size() - 1][0],100, tab, 10, F));
	vector<double> ans;
	for (int i = 0; i < (intervals.size() - 1) / 2; i = i+2)
	{
		ans.push_back(bisection(intervals[i], intervals[i+1], eps, tab, F));
	}
	
	return ans;
}


int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "«адача обратного интерполировани€." << endl;
	cout << "¬ведите число значений в таблице: ";
	int m;
	cin >> m;
	cout << endl << "¬ведите концы отрезка, из которого выбираютс€ узлы: ";
	double a, b;
	cin >> a >> b;
	cout << endl;
	vector<vector<double>> tab = table(a, b, m - 1);
	vector<vector<double>> rtab = reverse_table(a, b, m - 1);

	cout << "f(x)" << ' ' << "X(i)" << endl;
	for (int i = 0; i < rtab.size(); ++i)
	{
		cout << rtab[i][0] << ";  " << rtab[i][1] << endl;
	}
	double F;
	cout << endl << "¬ведите точку интерполировани€: ";
	cin >> F;
	int c;
	cout << endl << "¬ыберите способ обратного интерполировани€:";
	cin >> c;
	if (c == 1)
	{
		double x = backward_interpolation_w1(F, rtab, 10);
		cout << endl << "найденна€ точка:" << x << endl << "модуль нев€зки:" << abs(F - f(x)) << endl;
	}
	else if (c == 2)
	{
		vector<double> x2 = backward_interpolation_w2(F, tab, 10, 1e-6);
		cout << endl << "найденные точки:" << endl;
		for (int i = 0; i < x2.size(); ++i)
		{
			cout << x2[i] << ' ' << "модуль нев€зки:" << abs(F - f(x2[i])) << endl;
		}
		cout << endl;
	}
	return 0;
}