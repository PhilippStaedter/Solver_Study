#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_gardner1_Fig4B(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -C*X*reaction2_k1/pow(C + reaction2_K5, 2) + X*reaction2_k1/(C + reaction2_K5);
    dwdx[1] = reaction3_kd;
    dwdx[2] = -C*V1p*(-M + 1)/(pow(C + K6, 2)*(-M + reaction4_K1 + 1)) + V1p*(-M + 1)/((C + K6)*(-M + reaction4_K1 + 1));
    dwdx[3] = Y*reaction8_a1;
    dwdx[4] = C*V1p*(-M + 1)/((C + K6)*pow(-M + reaction4_K1 + 1, 2)) - C*V1p/((C + K6)*(-M + reaction4_K1 + 1));
    dwdx[5] = -M*reaction5_V2/pow(M + reaction5_K2, 2) + reaction5_V2/(M + reaction5_K2);
    dwdx[6] = V3p*(-X + 1)/(-X + reaction6_K3 + 1);
    dwdx[7] = C*reaction2_k1/(C + reaction2_K5);
    dwdx[8] = M*V3p*(-X + 1)/pow(-X + reaction6_K3 + 1, 2) - M*V3p/(-X + reaction6_K3 + 1);
    dwdx[9] = -X*reaction7_V4/pow(X + reaction7_K4, 2) + reaction7_V4/(X + reaction7_K4);
    dwdx[10] = reaction13_d1;
    dwdx[11] = C*reaction8_a1;
    dwdx[12] = reaction10_alpha*reaction10_d1;
    dwdx[13] = reaction11_alpha*reaction11_kd;
    dwdx[14] = reaction9_a2;
}