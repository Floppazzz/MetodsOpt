#include <cmath>
#include <iostream>
using namespace std;

double f(double x) { return max(x * x + 1.0, 2.0 * abs(x - 1.0)); }
