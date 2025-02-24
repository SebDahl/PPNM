#include <cmath>
#include <stdexcept>  // For std::invalid_argument

namespace sfuns{

double fgamma(double x) {
    if(x<0)return M_PI/sin(M_PI*x)/fgamma(1-x);
    if(x<9)return fgamma(x+1)/x;
    double lnfgamma=log(2*M_PI)/2+(x-0.5)*log(x)-x
        +(1.0/12)/x-(1.0/360)/(x*x*x)+(1.0/1260)/(x*x*x*x*x);
    return exp(lnfgamma);
    }

double lngamma(double x) {
    if (x <= 0) return NAN;  // Return NaN instead of throwing an exception
    if (x < 9) return lngamma(x + 1) - log(x);
    return x * log(x + 1 / (12 * x - 1 / (10 * x))) - x + log(2 * M_PI / x) / 2;
}
    

double erf(double x){
    /// single precision error function (Abramowitz and Stegun, from Wikipedia)
    double a1=0.278393, a2=0.230389, a3=0.000972, a4=0.078108;
    double y = fabs(x);
    double res = 1-1/pow((1+a1*y+a2*y*y+a3*y*y*y+a4*y*y*y*y),4);
    if (x<0) return -res;
    else return res;
}


}//class sfuns