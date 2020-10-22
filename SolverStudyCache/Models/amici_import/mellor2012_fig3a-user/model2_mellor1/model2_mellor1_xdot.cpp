#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model2_mellor1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -0.001*flux_r0 - 0.001*flux_r1 - 0.001*flux_r2;
    xdot[1] = 1.2e-5*flux_r0 + 1.5e-5*flux_r1 + 0.000107*flux_r2;
    xdot[2] = 0.00016200000000000001*flux_r0 + 0.000127*flux_r1 + 0.00021800000000000001*flux_r2;
    xdot[3] = 4.0000000000000003e-5*flux_r0 + 2.5999999999999998e-5*flux_r1 + 0.00021800000000000001*flux_r2;
    xdot[4] = 1.4e-5*flux_r0 + 1.8e-5*flux_r1 + 9.800000000000001e-5*flux_r2;
    xdot[5] = 3.9999999999999998e-6*flux_r0 + 1.5999999999999999e-5*flux_r1 + 9.7e-5*flux_r2;
    xdot[6] = 0.001*flux_r3 + 0.001*flux_r4;
    xdot[7] = 0.00057399999999999997*flux_r0 + 0.00075100000000000004*flux_r1 + 6.8000000000000013e-5*flux_r2 - 0.001*flux_r3;
    xdot[8] = 0.000144*flux_r0 + 2.3e-5*flux_r1 + 5.8999999999999998e-5*flux_r2 - 0.001*flux_r4;
    xdot[9] = 5.0000000000000002e-5*flux_r0 + 2.5000000000000001e-5*flux_r1 + 0.00013600000000000003*flux_r2;
}