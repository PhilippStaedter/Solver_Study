#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_borghans1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = -1.0*Rho*pow(Y, 2)*pow(Z, 8)*a*beta*v3_Vm3/(pow(d, 2)*(pow(Y, 2) + pow(v3_Ky, 2))*pow(pow(Z, 4)*a/d + 1, 2)) + 1.0*Rho*pow(Y, 2)*pow(Z, 4)*beta*v3_Vm3/(d*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
            break;
        case 1:
            dwdp[1] = 1.0*Rho*pow(Y, 2)*pow(Z, 4)*a*v3_Vm3/(d*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
            dwdp[6] = 1.0*vin_v1;
            break;
        case 2:
            dwdp[1] = 1.0*Rho*pow(Y, 2)*pow(Z, 8)*pow(a, 2)*beta*v3_Vm3/(pow(d, 3)*(pow(Y, 2) + pow(v3_Ky, 2))*pow(pow(Z, 4)*a/d + 1, 2)) - 1.0*Rho*pow(Y, 2)*pow(Z, 4)*a*beta*v3_Vm3/(pow(d, 2)*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
            break;
        case 3:
            dwdp[0] = -2.0*pow(Z, 2)*v2_K2*v2_Vm2/pow(pow(Z, 2) + pow(v2_K2, 2), 2);
            break;
        case 4:
            dwdp[0] = 1.0*pow(Z, 2)/(pow(Z, 2) + pow(v2_K2, 2));
            break;
        case 5:
            dwdp[1] = -2.0*Rho*pow(Y, 2)*pow(Z, 4)*a*beta*v3_Ky*v3_Vm3/(d*pow(pow(Y, 2) + pow(v3_Ky, 2), 2)*(pow(Z, 4)*a/d + 1));
            break;
        case 6:
            dwdp[1] = 1.0*Rho*pow(Y, 2)*pow(Z, 4)*a*beta/(d*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
            break;
        case 7:
            dwdp[2] = 1.0*Y;
            break;
        case 8:
            dwdp[3] = 1.0*Z;
            break;
        case 9:
            dwdp[4] = 1.0*Rho*pow(Z, 4);
            break;
        case 10:
            dwdp[5] = -1.0*Rho + 1.0;
            break;
        case 11:
            dwdp[6] = 1.0*beta;
            break;
        case 12:
            dwdp[6] = 1.0;
            break;
    }
}