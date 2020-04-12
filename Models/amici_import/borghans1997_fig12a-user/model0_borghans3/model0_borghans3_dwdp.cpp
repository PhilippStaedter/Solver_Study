#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_borghans3(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[3] = 1.0*pow(X, 2)*v3s_Vm3s/(pow(X, 2) + pow(v3s_K3s, 2));
            dwdp[7] = 1.0*vin_v1;
            break;
        case 1:
            dwdp[0] = -2.0*pow(Z, 2)*v2i_K2i*v2i_Vm2i/pow(pow(Z, 2) + pow(v2i_K2i, 2), 2);
            break;
        case 2:
            dwdp[0] = 1.0*pow(Z, 2)/(pow(Z, 2) + pow(v2i_K2i, 2));
            break;
        case 3:
            dwdp[1] = -2.0*pow(Z, 2)*v2s_K2s*v2s_Vm2s/pow(pow(Z, 2) + pow(v2s_K2s, 2), 2);
            break;
        case 4:
            dwdp[1] = 1.0*pow(Z, 2)/(pow(Z, 2) + pow(v2s_K2s, 2));
            break;
        case 5:
            dwdp[2] = -2.0*pow(Y, 2)*pow(Z, 2)*v3i_K3z*v3i_Vm3i/((pow(Y, 2) + pow(v3i_K3y, 2))*pow(pow(Z, 2) + pow(v3i_K3z, 2), 2));
            break;
        case 6:
            dwdp[2] = -2.0*pow(Y, 2)*pow(Z, 2)*v3i_K3y*v3i_Vm3i/(pow(pow(Y, 2) + pow(v3i_K3y, 2), 2)*(pow(Z, 2) + pow(v3i_K3z, 2)));
            break;
        case 7:
            dwdp[2] = 1.0*pow(Y, 2)*pow(Z, 2)/((pow(Y, 2) + pow(v3i_K3y, 2))*(pow(Z, 2) + pow(v3i_K3z, 2)));
            break;
        case 8:
            dwdp[3] = -2.0*pow(X, 2)*beta*v3s_K3s*v3s_Vm3s/pow(pow(X, 2) + pow(v3s_K3s, 2), 2);
            break;
        case 9:
            dwdp[3] = 1.0*pow(X, 2)*beta/(pow(X, 2) + pow(v3s_K3s, 2));
            break;
        case 10:
            dwdp[4] = 1.0*Y;
            break;
        case 11:
            dwdp[5] = 1.0*Z;
            break;
        case 12:
            dwdp[6] = 1.0*X;
            break;
        case 13:
            dwdp[7] = 1.0*beta;
            break;
        case 14:
            dwdp[7] = 1.0;
            break;
    }
}