#include <iostream>
#include <vector>
#include <cmath>
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

int main() {
    Point xk(1.0, 1.0);      
    double eps = 0.01;   
    
    // Параметры Армихо
    double alpha_init = 1.0;  
    double gamma = 0.1;       
    double theta = 0.5;       
    int k = 1;                
    int max_iter = 10000;     

    while (true) {
        Point grad = grad_f(xk);
        if (grad.norm() <= eps) {
            cout << "Сходимость достигнута на итерации " << k << endl;
            break;
        };

        Point pk = grad * (-1.0);
        double alpha_k = Armiho(alpha_init, gamma, theta, eps, xk, pk);
        xk = xk + pk * alpha_k;
        k++;
        if (k > max_iter) {
            cout << "Превышено максимальное число итераций" << endl;
            break;
        }
    }
    
    cout << "Точка минимума: (" << xk.x1 << ", " << xk.x2 << ")" << endl;
    cout << "Значение функции: " << f(xk) << endl;
    cout << "Норма градиента: " << grad_f(xk).norm() << endl;
    return 0;
}