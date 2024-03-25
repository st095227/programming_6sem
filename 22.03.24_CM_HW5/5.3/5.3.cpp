#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
double KF_Gauss(double f(double, double, vector<double>), double a, double b, double n, double k, vector<double> roots_node);

vector<double> Mom(double n, double c, double d)
{
    vector<double> ans;
    for (int i = 0; i < 2 * n; ++i)
    {
        ans.push_back(((double)2 / (3+2*i)) * (pow(d, 1.5+i) - pow(c, 1.5+i)));
    }
    return ans;
}
// Вывод системы уравнений
void sysout(vector<vector<double>> a, vector<double> y, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << a[i][j] << "*a" << j;
            if (j < n - 1)
                cout << " + ";
        }
        cout << " = " << y[i] << endl;
    }
    return;
}

vector<double> gauss(vector<vector<double>> a, vector<double> y, int n)
{
    vector<double> x;
    double max;
    int k, index;
    const double eps = 1e-12;  // точность
    for (int i = 0; i < n; ++i) x.push_back(0);
    k = 0;
    while (k < n)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(a[i][k]) > max)
            {
                max = abs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return { 0 };
        }
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = k; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    return x;
}

double poly(double x, vector<double> a)
{
    double sum = pow(x,a.size());
    for (int i = 0; i < a.size(); ++i)
    {
        sum += a[i] * pow(x, a.size() - i - 1);
    }
    return sum;
}

double secantPoly (double a, double b, double eps, vector<double> arr_ai)
{
    double x0 = a;
    double x1 = b;
    double xk = x1 - (poly(x1,arr_ai) / (poly(x1,arr_ai) - poly(x0,arr_ai))) * (x1 - x0);
    while (abs(xk - x1) > eps)
    {
        x0 = x1;
        x1 = xk;
        xk = x1 - (poly(x1, arr_ai) / (poly(x1, arr_ai) - poly(x0, arr_ai))) * (x1 - x0);
    }
    return xk;
}

vector<double> separate_roots(int N, double a, double b, vector<double> arr_ai)
{
    double H = (b-a) / N;
    double x1 = a;
    double x2 = x1 + H;
    double y1 = poly(x1, arr_ai);
    vector<double> arr;
    while (x2 <= 1)
    {
        double y2 = poly(x2, arr_ai);
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

vector<double> node_roots(double a, double b, vector<double> arr_ai)
{
    vector<double> segments;
    vector<double> node_r;
    segments = separate_roots(1000, a, b, arr_ai);
    for (int i = 0; i < segments.size(); i = i + 2)
    {
        node_r.push_back(secantPoly(segments[i], segments[i + 1], 1e-12, arr_ai));
    }
    return node_r;
}

double func_coef(double x, double k, vector<double> roots)
{
    double ans = sqrt(x);
    for (int i = 0; i < roots.size(); ++i)
    {
        if (i != k)
        {
            ans *= (x - roots[i]);
        }
    }
    for (int i = 0; i < roots.size(); ++i)
    {
        if (i != k)
        {
            ans /= (roots[k] - roots[i]);
        }
    }
    return ans;
}

vector<double> coef(vector<double> roots, double a, double b)
{
    vector<double> ans;
    for (double i = 0; i < roots.size(); ++i)
    {
        ans.push_back(KF_Gauss(func_coef,a,b,roots.size(),i,roots));
    }

    return ans;
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
        return (((2 * n - 1) / n) * Legendre(n - 1, x) * x - (((n - 1) / n) * Legendre(n - 2, x)));
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

vector<double> separate_rootsLeg(int N, int n)
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

double KF_Gauss(double f(double, double, vector<double>), double a, double b, double n, double k, vector<double> roots_node)
{
    vector<double> segm;
    vector<double> roots;
    vector<double> coeff;
    segm = separate_rootsLeg(1000, n);
    for (int i = 0; i < segm.size(); i = i + 2)
    {
        roots.push_back(secantLeg(segm[i], segm[i + 1], 1e-12, n));
    }
    for (int i = 0; i < roots.size(); ++i)
    {
        coeff.push_back(((2 * (1 - roots[i] * roots[i])) / (n * n * Legendre(n - 1, roots[i]) * Legendre(n - 1, roots[i]))));
    }
    double sum = 0;
    for (int i = 0; i < coeff.size(); ++i)
    {
        sum += coeff[i] * f(0.5 * (a + b + (b - a) * roots[i]), k, roots_node);
    }
    sum *= (0.5 * (b - a));
    return sum;
}

double KF_Gauss_Type(double a, double b, double N, double m)
{
    double h = (b - a) / m;
    double ans = 0;
    double sum = 0;
    vector<double> mom;
    vector<vector<double>> sys;
    vector<double> string;
    vector<double> y;
    for (int i = 1; i <= m; ++i)
    {
        mom = Mom(N, a + (i - 1) * h, a + i * h);
        for (int j = N; j < 2 * N; ++j)
        {
            y.push_back(-mom[j]);
        }
        for (int j = N-1; j < (2 * N) - 1; ++j)
        {
            for (int k = j; k > j-N; --k)
            {
                string.push_back(mom[k]);
            }
            sys.push_back(string);
            string.clear();
        }
        vector<double> ai = gauss(sys, y, N);
        vector<double> node = node_roots(a + (i - 1) * h, a + i * h, ai);
        vector<double> coeff = coef(node, a + (i - 1) * h, a + i * h);
        for (int k = 0; k < N; ++k)
        {
            sum += coeff[k] * sin(node[k]);
        }
        ans += sum;
        sum = 0;
        y.clear();
        sys.clear();
    }
    return ans;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Приближённое вычисление интеграла при помощи составной КФ Гаусса." << endl;
    cout << "Введите промежуток интегрирования: ";
    double a, b;
    cin >> a >> b;
    cout << endl << "Введите количество узлов: ";
    int N;
    cin >> N;
    cout << endl << "Введите количество разбиений: ";
    int m;
    cin >> m;
    cout << KF_Gauss_Type(a, b, N, m);
    return 0;
}