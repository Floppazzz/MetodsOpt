#include <cmath>
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

double Fibonachi(double eps, double ak, double bk, double alpha = 1e-6) {
	vector<double> F = {0, 1};
    int k = 1;
    
    while ((bk - ak) / F[k] > eps) {
        F.push_back(F[k] + F[k - 1]);
        k++;
    }
    int n = k;
    
    double xk = ak + (F[n - 2] / F[n]) * (bk - ak);
    double yk = ak + (F[n - 1] / F[n]) * (bk - ak);
    double fx = f(xk);
    double fy = f(yk);

    for (int i = 1; i < n - 1; i++) {
        if (fx > fy) {
            ak = xk;  
            xk = yk;  
            fx = fy;
            yk = ak + (F[n - i - 1] / F[n - i]) * (bk - ak);
            fy = f(yk);
        } else {
            bk = yk;  
            yk = xk;  
            fy = fx;
            xk = ak + (F[n - i - 2] / F[n - i]) * (bk - ak);
            fx = f(xk);
        }
    }

    double xn = xk;
    double yn = yk + alpha;
    
    if (f(xn) > f(yn)) {
        ak = xn;
    } else {
        bk = yn;
    }

    return (ak + bk) / 2.0;
}


int main(){
    double eps = 0.01;
    double delta = (sqrt(5.0) - 1.0) / 2.0;
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
