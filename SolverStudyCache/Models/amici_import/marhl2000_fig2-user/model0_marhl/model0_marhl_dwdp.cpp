#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_marhl(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -2.0*pow(Ca_cyt, 2)*v1_K1*v1_Kch*(CaER - Ca_cyt)/pow(pow(Ca_cyt, 2) + pow(v1_K1, 2), 2);
            break;
        case 1:
            dwdp[0] = 1.0*pow(Ca_cyt, 2)*(CaER - Ca_cyt)/(pow(Ca_cyt, 2) + pow(v1_K1, 2));
            break;
        case 2:
            dwdp[1] = 1.0*CaPr;
            break;
        case 3:
            dwdp[2] = 1.0*Ca_cyt*Pr;
            break;
        case 4:
            dwdp[3] = 1.0*CaER - 1.0*Ca_cyt;
            break;
        case 5:
            dwdp[4] = 1.0*Ca_cyt;
            break;
        case 6:
            dwdp[5] = 1.0*CaM;
            break;
        case 7:
            dwdp[5] = -2.0*CaM*pow(Ca_cyt, 2)*v7_K3*v7_Kout/pow(pow(Ca_cyt, 2) + pow(v7_K3, 2), 2);
            break;
        case 8:
            dwdp[5] = 1.0*CaM*pow(Ca_cyt, 2)/(pow(Ca_cyt, 2) + pow(v7_K3, 2));
            break;
        case 9:
            dwdp[6] = -8.0*pow(Ca_cyt, 8)*pow(v9_K2, 7)*v9_Kin/pow(pow(Ca_cyt, 8) + pow(v9_K2, 8), 2);
            break;
        case 10:
            dwdp[6] = 1.0*pow(Ca_cyt, 8)/(pow(Ca_cyt, 8) + pow(v9_K2, 8));
            break;
    }
}