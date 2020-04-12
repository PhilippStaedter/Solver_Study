#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_mellor1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -0.001*dwdx0 - 0.001*dwdx1 - 0.001*dwdx2;
    J[10] = 1.2e-5*dwdx0 + 1.5e-5*dwdx1 + 0.000107*dwdx2;
    J[20] = 0.00016200000000000001*dwdx0 + 0.000127*dwdx1 + 0.00021800000000000001*dwdx2;
    J[30] = 4.0000000000000003e-5*dwdx0 + 2.5999999999999998e-5*dwdx1 + 0.00021800000000000001*dwdx2;
    J[40] = 1.4e-5*dwdx0 + 1.8e-5*dwdx1 + 9.800000000000001e-5*dwdx2;
    J[50] = 3.9999999999999998e-6*dwdx0 + 1.5999999999999999e-5*dwdx1 + 9.7e-5*dwdx2;
    J[67] = 0.001*dwdx3;
    J[68] = 0.001*dwdx4;
    J[70] = 0.00057399999999999997*dwdx0 + 0.00075100000000000004*dwdx1 + 6.8000000000000013e-5*dwdx2;
    J[77] = -0.001*dwdx3;
    J[80] = 0.000144*dwdx0 + 2.3e-5*dwdx1 + 5.8999999999999998e-5*dwdx2;
    J[88] = -0.001*dwdx4;
    J[90] = 5.0000000000000002e-5*dwdx0 + 2.5000000000000001e-5*dwdx1 + 0.00013600000000000003*dwdx2;
}