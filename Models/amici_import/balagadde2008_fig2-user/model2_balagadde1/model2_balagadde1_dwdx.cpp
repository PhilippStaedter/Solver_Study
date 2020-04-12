#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_balagadde1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*C2*(-2*pow(A1, 3)*d2/pow(pow(A1, 2) + K2, 2) + 2*A1*d2/(pow(A1, 2) + K2));
    dwdx[1] = 1.0*D + 1.0*dAA1;
    dwdx[2] = -2.0*A2*C1*K1*(pow(IPTG, 2)/(pow(IPTG, 2) + 25) + 0.5)/pow(pow(A2, 2) + K1, 2);
    dwdx[3] = 1.0*D + 1.0*dAA2;
    dwdx[4] = -1.0*C1*kc1/Cm + 1.0*kc1*(1 - (C1 + C2)/Cm);
    dwdx[5] = 1.0*D + 1.0*K1*(pow(IPTG, 2)/(pow(IPTG, 2) + 25) + 0.5)/(pow(A2, 2) + K1);
    dwdx[6] = -1.0*C2*kc2/Cm;
    dwdx[7] = 1.0*kA1;
    dwdx[8] = -1.0*C1*kc1/Cm;
    dwdx[9] = -1.0*C2*kc2/Cm + 1.0*kc2*(1 - (C1 + C2)/Cm);
    dwdx[10] = 1.0*pow(A1, 2)*d2/(pow(A1, 2) + K2) + 1.0*D;
    dwdx[11] = 0.029999999999999999*pow(IPTG, 2)/(pow(IPTG, 2) + 25) + 0.02;
    dwdx[12] = 1.0*C1*K1*(-2*pow(IPTG, 3)/pow(pow(IPTG, 2) + 25, 2) + 2*IPTG/(pow(IPTG, 2) + 25))/(pow(A2, 2) + K1);
    dwdx[13] = 1.0*C2*(-0.059999999999999998*pow(IPTG, 3)/pow(pow(IPTG, 2) + 25, 2) + 0.059999999999999998*IPTG/(pow(IPTG, 2) + 25));
}