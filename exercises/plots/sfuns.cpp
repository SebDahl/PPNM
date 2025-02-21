#include <cmath>
namespace sfuns{

double fgamma(double x){
///single precision gamma function (formula from Wikipedia)
if(x<0)return M_PI/sin(M_PI*x)/fgamma(1-x); // Euler's reflection formula
if(x<9)return fgamma(x+1)/x; // Recurrence relation
double lnfgamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
return exp(lnfgamma);
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