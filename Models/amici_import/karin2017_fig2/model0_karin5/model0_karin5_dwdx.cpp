#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_karin5(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 0.030000000000000002/(pow(alpha/g, 1.7) + 1);
    dwdx[1] = -1.0/1440.0*tamox + 6.9444444444444444e-5/(37.262740006283494*pow(1.0/g, 1.7) + 1) - 0.00013888888888888889/((1.0/65536.0)*pow(g, 8) + 1) - 0.00013888888888888889/(1 + 11390625/pow(g, 6));
    dwdx[2] = (1.0/1440.0)*tamox;
    dwdx[3] = -0.00050000000000000001*ins - 0.001;
    dwdx[4] = 0.051000000000000004*beta*pow(alpha/g, 1.7)/(g*pow(pow(alpha/g, 1.7) + 1, 2)) + 0.051000000000000004*mbeta*pow(alpha/(amici_k*g), 1.7)/(g*pow(pow(alpha/(amici_k*g), 1.7) + 1, 2));
    dwdx[5] = (1.0/1440.0)*beta*(2.4414062500000001e-5*pow(g, 7)/pow((1.0/65536.0)*pow(g, 8) + 1, 2) + 6.3346658010681942*pow(1.0/g, 1.7)/(g*pow(37.262740006283494*pow(1.0/g, 1.7) + 1, 2)) - 13668750.0/(pow(g, 7)*pow(1 + 11390625/pow(g, 6), 2)));
    dwdx[6] = (1.0/1440.0)*mbeta*(2.4414062500000001e-5*pow(amici_k, 8)*pow(g, 7)/pow((1.0/65536.0)*pow(amici_k, 8)*pow(g, 8) + 1, 2) + 6.3346658010681942*pow(1/(amici_k*g), 1.7)/(g*pow(37.262740006283494*pow(1/(amici_k*g), 1.7) + 1, 2)) - 13668750.0/(pow(amici_k, 6)*pow(g, 7)*pow(1 + 11390625/(pow(amici_k, 6)*pow(g, 6)), 2)));
    dwdx[7] = -0.00050000000000000001*g;
    dwdx[8] = -3.0/10.0;
    dwdx[9] = 0.030000000000000002/(pow(alpha/(amici_k*g), 1.7) + 1);
    dwdx[10] = -0.00013888888888888889/((1.0/65536.0)*pow(amici_k, 8)*pow(g, 8) + 1) + 6.9444444444444444e-5/(37.262740006283494*pow(1/(amici_k*g), 1.7) + 1) - 0.00013888888888888889/(1 + 11390625/(pow(amici_k, 6)*pow(g, 6)));
    dwdx[11] = -1.0/1440.0*beta;
    dwdx[12] = (1.0/1440.0)*beta;
    dwdx[13] = -0.0010416666666666667*M_LN2;
}