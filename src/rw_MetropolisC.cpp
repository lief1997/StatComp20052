#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
float lap_f(int x) {
  float lap = exp(-abs(x));
  return lap;
}
// [[Rcpp::export]]
List rw_MetropolisC(double sigma, double x0, int N) {
  NumericVector x(N);
  NumericVector u(N);
  x[0] = x0;
  u = runif(N);
  int k = 0;
  for (int i = 1; i < N; i++) {
    double y;
    y = rnorm(1,x[i-1],sigma)[0];
    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) 
      x[i] = y;
    else {
      x[i] = x[i-1];
      k++;
    }
  }
  List out;
  out["x"]=x;
  out["k"]=k;
  return out;
  //return x;
}

