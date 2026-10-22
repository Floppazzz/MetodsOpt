#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
double f(double x){
	return (fabs(x-3)+fabs(x+2)+fabs(x+1));
}
double method1(const double&a_,const double&b_,const double&eps,std::function<double(double)>f){
	double a=a_;
	double b=b_;
	double t = (sqrt(5.0)-1)/2;
	double x=a+(1-t)*(b-a);
	double y=a+t*(b-a);
	double fx = f(x);
	double fy=f(y);
	
	while ((b-a)>eps){
		std::cout<<x<<" "<<y<<" "<<f(x)<<" "<<f(y)<<std::endl;
		if(fx>fy){
			a=x;
			//b=b;
			x=y;
			fx=fy;
			y=a+t*(b-a);
			fy=f(y);
		}
		else{
			//a=a;
			b=y;
			y=x;
			fy=fx;
			x=a+(1-t)*(b-a);
			fx=f(x);
		}
	} return (a+b)/2;
}
double method2(const double& a_,const double& b_, const double& eps, std::function<double(double)>f){
    double a=a_;double b=b_;double alpha=eps/3;
	std::vector<double> F = {1, 1};
    int k = 1;
    
    while ((b - a) / F[k] > eps) {
        F.push_back(F[k] + F[k - 1]);
        k++;
    }
    int n = k;

    double x = a + (F[n - 2] / F[n]) * (b - a);
    double y = a + (F[n - 1] / F[n]) * (b - a);
    double fx = f(x);
    double fy = f(y);

    for (int i = 1; i < n - 1; i++) {
        if (fx > fy) {
            a = x;  
            x = y;  
            fx = fy;
            y = a + (F[n - i - 1] / F[n - i]) * (b - a);
            fy = f(y);
        } else {
            b = y;  
            y = x;  
            fy = fx;
            x = a + (F[n - i - 2] / F[n - i]) * (b - a);
            fx = f(x);
        }
    }

    double xn = x;
    double yn = y + alpha;
    
    if (f(xn) > f(yn)) {
        a = xn;
    } else {
        b = yn;
    }

    return (a + b) / 2.0;
}

int main(int argc, char** argv){
	std::cout<<method1(-3,4,0.01,f)<<"\n"<<method2(-3,4,0.01,f)<<std::endl;
}