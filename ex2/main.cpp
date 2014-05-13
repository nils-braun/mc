// Exercise 2 for MC-Ereignisgeneratoren
// May 2014
//

#include <stdlib.h>
#include <complex>
#include <cmath>

const double G_F = 10e-5;
const double M_W = 80;
const double Gamma_W = 2;

double calcSumM(const double Theta, const double sHat) {
  return((2*G_F*pow(M_W, 8.0) * pow(1 + cos(Theta), 2.0))/(pow(sHat - pow(M_W, 2), 2.0) + pow(M_W, 2.0)*pow(Gamma_W, 2.0)));
}

double pdf(const double x, const double sHat) {
  return(pow(x, 0.2-0.3*log(sHat)) * pow(1-x, 0.1));
}

int main() {

}
