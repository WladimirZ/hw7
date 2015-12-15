#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void func (double* const y0){
  
  double y[4] = { y0[0], y0[1], y0[2], y0[3] };
  
  const double mu = 0.012277471;
  double r = sqrt((y[0]+mu)*(y[0]+mu)+y[2]*y[2]);
  double s = sqrt((y[0]-1+mu)*(y[0]-1+mu)+y[2]*y[2]);
  
  y0[0] = y[1];
  y0[1] = y[0] + 2.0*y[3] - (1.0 - mu) * (y[0] + mu) / (r*r*r) - mu * (y[0] - 1 + mu) / (s*s*s);
  y0[2] = y[3];
  y0[3] = y[2] - 2.0*y[1] - (1.0 - mu) * y[2] / (r*r*r) - mu*y[2]/(s*s*s);
  
}

void RK (double* const y4, double* const y5, const double dx) { // dx = dt
  const int d = 4;
  double k1[d], k2[d], k3[d], k4[d], k5[d], k6[d], k7[d];
  
  for (int i = 0; i<d; i++) k1[i] = y4[i];
  func(k1);
  
  for (int i = 0; i<d; i++) k2[i] = y4[i] + 0.2 * dx * k1[i];
  func(k2);
  
  for (int i = 0; i<d; i++) k3[i] = y4[i] + 3.0/40 * dx * k1[i] + 9.0/40 * dx * k2[i];
  func(k3);
  
  for (int i = 0; i<d; i++) k4[i] = y4[i] + 44.0/45 * dx * k1[i] - 56.0/15 * dx * k2[i] + 32.0/9 * dx *k3[i];
  func(k4);
  
  for (int i = 0; i<d; i++) k5[i] = y4[i] + 19372.0/6561 * dx * k1[i] - 25360.0/2187 * dx * k2[i] + 64448.0/6561 * dx * k3[i] - 212.0/729 * dx * k4[i];
  func(k5);
  
  for (int i = 0; i<d; i++) k6[i] = y4[i] + 9017.0/3168 * dx * k1[i] - 355.0/33 *dx * k2[i] + 46732.0/5247 * dx * k3[i] + 49.0/176 * dx * k4[i] - 5103.0/18656 * dx * k5[i];
  func(k6);
  
  for (int i = 0; i<d; i++) k7[i] = y4[i] + 35.0/384 * dx * k1[i] + 500.0/1113 * dx * k3[i] + 125.0/192 * dx * k4[i] - 2187.0/6784 * dx * k5[i] + 11.0/84 * dx * k6[i];
  func(k7);
  
  for (int i = 0; i<d; i++) y4[i] = y4[i] + 5179.0/57600 * dx * k1[i] + 7571.0/16695 * dx * k3[i] + 393.0/640 * dx * k4[i] - 92097.0/339200 * dx *k5[i] + 187.0/2100 * dx * k6[i] + 1.0/40 * dx * k7[i];
  
  for (int i = 0; i<d; i++) y5[i] = y4[i] + 35.0/384 * dx * k1[i] + 500.0/1113 * dx * k3[i] + 125.0/192 * dx * k4[i] - 2187.0/6784 * dx * k5[i] + 11.0/84 * dx * k6[i];
  
}

double error ( double* y4, double* y5) {
  return fmax ( fmax (abs(y4[0]-y5[0]), abs(y4[1]-y5[1]) ), fmax( abs(y4[2]-y5[2]), abs(y4[3]-y5[3]) ) );
}

int main() {
  double dt = 1e-5;
  const double Tol = 1e-5;
  double Xi;
  double y4[4] = {0.994, 0, 0, -2.00158510637908};
  double y5[4] = {0.994, 0, 0, -2.00158510637908};
  
  ofstream out ( "sol" );
  double t;
  for ( t = 0.0; t<30; t += dt) {
    out << t << "\t" << dt << "\t" << y4[0]  << "\t" << y4[2] << endl;
    RK ( y4, y5, dt );
    Xi = error (y4, y5);
    dt *= 0.9 * pow( Tol/Xi, 1.0/5 );
  }
  out << t << "\t" << dt << "\t" << y4[0]  << "\t" << y4[2] << endl;
  out.close();
  
  return 0;
}