#include <iostream>
#include <cmath>
#include <functional>
using namespace std;

class Point {
    public:
        double x1, x2;
        Point(double x1, double x2) : x1(x1), x2(x2) {}

        Point operator+(const Point& other) const {
            return Point(x1 + other.x1, x2 + other.x2);
        }
        Point operator*(double scalar) const {
            return Point(x1 * scalar, x2 * scalar);
        }
        double norm() const {
            return sqrt(x1 * x1 + x2 * x2);
        }
};

Point operator*(double scalar, const Point& point) {
    return point * scalar;
}

double f(const Point& p) {
    // f(x1, x2) = exp(x1) + (x1 + x2)^4 + exp(2x2)
    return exp(p.x1) + pow(p.x1 + p.x2, 4) + exp(2 * p.x2);
}

Point grad_f(const Point& p) {
    const double df_dx1 = exp(p.x1) + 4 * pow(p.x1 + p.x2, 3);
    const double df_dx2 = 4 * pow(p.x1 + p.x2, 3) + 2 * exp(2 * p.x2);
    return Point(df_dx1, df_dx2);
}

double Armiho(double alpha, double gamma, double theta, double eps, const Point& xk, const Point& pk) {
    double fxk = f(xk);
    Point grad = grad_f(xk);
    double grad_norm_sqr = grad.norm() * grad.norm();

    while (true) {
        Point xk_new = xk + pk * alpha;
        double fxk_new = f(xk_new);
        
        if (fxk_new - fxk <= -gamma * alpha * grad_norm_sqr) {
            return alpha;
        }
        alpha = theta * alpha;
    }
}

double gold_section(function<double(double)> func, double ak, double bk, double eps, double delta){
    double xk = ak + (1 - delta)*(bk - ak);
    double yk = ak + delta * (bk-ak);
    double fx = func(xk);
    double fy = func(yk);
    while (bk - ak > eps){
        if (fx > fy){
            ak = xk;
            xk = yk;
            fx = fy;
            yk = ak + delta * (bk - ak);
            fy = func(yk);
        }
        else {
            bk = yk;
            yk = xk;
            fy = fx;
            xk = ak + (1-delta)*(bk -ak);
            fx = func(xk);
        }
    }
    return (ak + bk) / 2.0;
}


int main() {
    Point xk_armiho(1.0, 1.0);
    Point xk_gold(1.0, 1.0);      
    double eps = 0.01;   
    double delta = (sqrt(5.0) - 1.0) / 2.0;
    double alpha_max = 5.0;
    
    // Параметры Армихо
    double alpha_init = 1.0;  
    double gamma = 0.1;       
    double theta = 0.5;       
    int k_armiho = 1;
    int k_gold = 1;                
    int max_iter = 10000;     

    while (true) {
        Point grad = grad_f(xk_armiho);
        if (grad.norm() <= eps) {
            cout << "Сходимость достигнута на итерации " << k_armiho << endl;
            break;
        };

        Point pk = grad * (-1.0);
        double alpha_k = Armiho(alpha_init, gamma, theta, eps, xk_armiho, pk);
        xk_armiho = xk_armiho + pk * alpha_k;
        k_armiho++;
        if (k_armiho > max_iter) {
            cout << "Превышено максимальное число итераций" << endl;
            break;
        }
    }

    while (true) {
        Point grad = grad_f(xk_gold);
        if (grad.norm() <= eps) {
            cout << "Сходимость достигнута на итерации " << k_gold << endl;
            break;
        }
        
        auto F = [xk_gold, grad] (double alpha) {
            Point next_x = xk_gold + grad * (-alpha);
            return f(next_x);
        };

        double alpha_k = gold_section(F, 0.0, alpha_max, eps, delta);
        xk_gold = xk_gold + grad * (-alpha_k);
        k_gold++;
        if (k_gold > max_iter) {
            cout << "Превышено максимальное число итераций" << endl;
            break;
        }
    }
    cout << "Метод градиентного спуска с использованием правила Армихо" << endl;
    cout << "Точка минимума: (" << xk_armiho.x1 << ", " << xk_armiho.x2 << ")" << endl;
    cout << "Значение функции: " << f(xk_armiho) << endl;
    cout << "Норма градиента: " << grad_f(xk_armiho).norm() << endl;
    cout << "Метод градиентного спуска с использованием золотого сечения" << endl;
    cout << "Количество итераций: " << k_gold << endl;
    cout << "Точка минимума: (" << xk_gold.x1 << ", " << xk_gold.x2 << ")" << endl;
    cout << "Значение функции: " << f(xk_gold) << endl;
    cout << "Норма градиента: " << grad_f(xk_gold).norm() << endl;
    return 0;
}