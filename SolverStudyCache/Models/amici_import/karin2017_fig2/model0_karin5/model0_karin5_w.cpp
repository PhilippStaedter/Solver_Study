#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_karin5(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = -g*(0.00050000000000000001*ins + 0.001) + 1.0/30.0;
    w[1] = 0.030000000000000002*beta/(pow(alpha/g, 1.7) + 1) - 3.0/10.0*ins + 0.030000000000000002*mbeta/(pow(alpha/(amici_k*g), 1.7) + 1);
    w[2] = (1.0/1440.0)*beta*(-tamox + 0.10000000000000001/(37.262740006283494*pow(1.0/g, 1.7) + 1) - 0.20000000000000001/((1.0/65536.0)*pow(g, 8) + 1) - 0.20000000000000001/(1 + 11390625/pow(g, 6)));
    w[3] = (1.0/1440.0)*beta*tamox + (1.0/1440.0)*mbeta*(-0.20000000000000001/((1.0/65536.0)*pow(amici_k, 8)*pow(g, 8) + 1) + 0.10000000000000001/(37.262740006283494*pow(1/(amici_k*g), 1.7) + 1) - 0.20000000000000001/(1 + 11390625/(pow(amici_k, 6)*pow(g, 6))));
    w[4] = -0.0010416666666666667*tamox*M_LN2;
}