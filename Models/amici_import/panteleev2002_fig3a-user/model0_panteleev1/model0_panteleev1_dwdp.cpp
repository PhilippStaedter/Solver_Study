#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_panteleev1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*VIIa_TF*X;
            break;
        case 1:
            dwdp[0] = -1.0*VIIa_TF_X;
            break;
        case 2:
            dwdp[1] = 1.0*VIIa_TF_X;
            break;
        case 3:
            dwdp[2] = 1.0*VIIa_TF_Xa;
            break;
        case 4:
            dwdp[2] = -1.0*VIIa_TF*Xa;
            break;
        case 5:
            dwdp[3] = 1.0*TFPI*Xa;
            break;
        case 6:
            dwdp[3] = -1.0*Xa_TFPI;
            break;
        case 7:
            dwdp[4] = 1.0*VIIa_TF*Xa_TFPI;
            break;
        case 8:
            dwdp[4] = -1.0*Xa_TFPI_VIIa_TF;
            break;
        case 9:
            dwdp[5] = 1.0*TFPI*VIIa_TF_Xa;
            break;
        case 10:
            dwdp[5] = -1.0*VIIa_TF_Xa_TFPI;
            break;
        case 11:
            dwdp[6] = 1.0*VIIa_TF_X*Xa_TFPI;
            break;
        case 12:
            dwdp[6] = -1.0*VIIa_TF_Xa_TFPI*X;
            break;
        case 13:
            dwdp[7] = 1.0*VIIa_TF_Xa_TFPI;
            break;
        case 14:
            dwdp[7] = -1.0*Xa_TFPI_VIIa_TF;
            break;
    }
}