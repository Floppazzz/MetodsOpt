#include <cmath>
#include <iostream>
using namespace std;

double f(double x) { return max(x * x + 1.0, 2.0 * abs(x - 1.0)); }

double dichotomy(double a, double b, double eps, double delta) {
  double ak = a;
  double bk = b;
  int k = 1;
  while (bk - ak > eps) {
    double x1 = (ak + bk) / 2.0 - delta;
    double x2 = (ak + bk) / 2.0 + delta;
    if (f(x1) > f(x2)) {
      ak = x1;
    } else {
      bk = x2;
    }
    k++;
  }
  return (ak + bk) / 2.0;
}

double halfmetod(double a, double b, double eps) {
  double ak = a;
  double bk = b;
  int k = 1;
  double ck = (ak + bk) / 2.0;
  double fc = f(ck);
  while (bk - ak > eps) {
    double x1 = ak + (bk - ak) / 4.0;
    double x2 = bk - (bk - ak) / 4.0;
    if (f(x1) < fc) {
      bk = ck;
      ck = x1;
      fc = f(x1);
    } else if (f(x2) < fc) {
      ak = ck;
      ck = x2;
      fc = f(x2);
    } else {
      ak = x1;
      bk = x2;
    }
    k++;
  }
  return ck;
}

int main() {
  double a = -4.0;
  double b = 3.0;
  double eps = 0.001;
  double delta = 0.0001;
  double first_metod = dichotomy(a, b, eps, delta);
  double second_metod = halfmetod(a, b, eps);
  cout << "Метод дихотомии: " << first_metod << endl;
  cout << "Метод половинного деления: " << second_metod << endl;
  return 0;
}
