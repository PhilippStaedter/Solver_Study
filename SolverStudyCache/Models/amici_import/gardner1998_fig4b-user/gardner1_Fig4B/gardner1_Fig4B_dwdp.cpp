#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_gardner1_Fig4B(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[7] = -C*V1p*(-M + 1)/(pow(C + K6, 2)*(-M + reaction4_K1 + 1));
            break;
        case 1:
            dwdp[7] = C*(-M + 1)/((C + K6)*(-M + reaction4_K1 + 1));
            break;
        case 2:
            dwdp[9] = M*(-X + 1)/(-X + reaction6_K3 + 1);
            break;
        case 3:
            dwdp[0] = 1;
            break;
        case 4:
            dwdp[1] = Z*reaction10_alpha;
            break;
        case 5:
            dwdp[1] = Z*reaction10_d1;
            break;
        case 6:
            dwdp[2] = Z*reaction11_kd;
            break;
        case 7:
            dwdp[2] = Z*reaction11_alpha;
            break;
        case 8:
            dwdp[3] = 1;
            break;
        case 9:
            dwdp[4] = Y;
            break;
        case 10:
            dwdp[5] = -C*X*reaction2_k1/pow(C + reaction2_K5, 2);
            break;
        case 11:
            dwdp[5] = C*X/(C + reaction2_K5);
            break;
        case 12:
            dwdp[6] = C;
            break;
        case 13:
            dwdp[7] = -C*V1p*(-M + 1)/((C + K6)*pow(-M + reaction4_K1 + 1, 2));
            break;
        case 14:
            dwdp[8] = -M*reaction5_V2/pow(M + reaction5_K2, 2);
            break;
        case 15:
            dwdp[8] = M/(M + reaction5_K2);
            break;
        case 16:
            dwdp[9] = -M*V3p*(-X + 1)/pow(-X + reaction6_K3 + 1, 2);
            break;
        case 17:
            dwdp[10] = X/(X + reaction7_K4);
            break;
        case 18:
            dwdp[10] = -X*reaction7_V4/pow(X + reaction7_K4, 2);
            break;
        case 19:
            dwdp[11] = C*Y;
            break;
        case 20:
            dwdp[12] = Z;
            break;
    }
}