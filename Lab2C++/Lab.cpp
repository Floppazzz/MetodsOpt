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

double Fibonachi(double eps, double a, double b, double alpha = 1e-6) {
    if (b - a <= eps) {
        return (a + b) / 2.0;
    }

    vector<double> F = {0.0, 1.0};
    while (F.back() < 2.0 * (b - a) / eps) {
        F.push_back(F[F.size() - 1] + F[F.size() - 2]);
    }

    const int n = static_cast<int>(F.size()) - 1;
    

    double xk = a + F[n - 2] / F[n] * (b - a);
    double yk = a + F[n - 1] / F[n] * (b - a);
    
    for (int k = n; k >= 3; --k) {
        double fx = f(xk);
        double fy = f(yk);
        
        if (fx > fy) {
            a = xk;
            xk = yk;
            yk = a + F[k - 2] / F[k - 1] * (b - a);
        } else {
            b = yk;
            yk = xk;
            xk = a + F[k - 3] / F[k - 1] * (b - a);
        }
    }
    
    double xn = xk;
    double yn = yk;
    
    if (abs(xn - yn) < alpha) {
        yn = xn + alpha;
    }
    
    double fxn = f(xn);
    double fyn = f(yn);
    
    if (fxn > fyn) {
        a = xn;
    } else {
        b = yn;
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

    const double xminF = Fibonachi(eps, ak, bk);
    const double xminGold = gold_section(ak, bk, eps, delta);

    cout << "Метод золотого сечения:" << endl;
    cout << "Точка минимума: " << xminGold << endl;
    cout << "Значение функции: " << f(xminGold) << endl;

    cout << "Метод Фибоначчи:" << endl;
    cout << "Точка минимума: " << xminF << endl;
    cout << "Значение функции: " << f(xminF) << endl;
}
