#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

double f(double x) { return max(x * x + 1.0, 2.0 * abs(x - 1.0)); }

double gold_section(double ak, double bk, double eps, double delta){
    double xk = ak + (1- delta)*(bk - ak);
    double yk = ak + delta * (bk-ak);
    double fx = f(xk);
    double fy = f(yk);
    while (bk - ak > eps){
        if (fx > fy){
            ak = xk;
            xk = yk;
            fx = fy;
            yk = ak + delta * (bk - ak);
            fy = f(yk);
        }
        else {
            bk = yk;
            yk = xk;
            fy = fx;
            xk = ak + (1-delta)*(bk -ak);
            fx = f(xk);
        }
    }
    return (ak + bk) / 2.0;
}


double Fibonachi(double eps, double a, double b){
    if (b - a <= eps) {
        return (a + b) / 2.0;
    }

    vector<double> fib = {0.0, 1.0};
    while (fib.back() < 2.0 * (b - a) / eps) {
        fib.push_back(fib[fib.size() - 1] + fib[fib.size() - 2]);
    }

    const int n = static_cast<int>(fib.size()) - 1;
    for (int k = n; k >= 4; --k) {
        double xk = a + fib[k - 2] / fib[k] * (b - a);
        double yk = a + fib[k - 1] / fib[k] * (b - a);
        double fx = f(xk);
        double fy = f(yk);

        if (fx > fy) {
            a = xk;
        } else {
            b = yk;
        }
    }

    return (a + b) / 2.0;
}

int main(){
    const double eps = 0.01;
    const double delta = (sqrt(5.0) - 1.0) / 2.0;
    double ak;
    double bk;

    cout << "Введите границы a и b:" << endl;
    cin >> ak >> bk;

    const double xminFib = Fibonachi(eps, ak, bk);
    const double xminGold = gold_section(ak, bk, eps, delta);

    cout << "Метод золотого сечения:" << endl;
    cout << "Точка минимума: " << xminGold << endl;
    cout << "Значение функции: " << f(xminGold) << endl;

    cout << "Метод Фибоначчи:" << endl;
    cout << "Точка минимума: " << xminFib << endl;
    cout << "Значение функции: " << f(xminFib) << endl;
}