#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_gorlich1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.2000000000000001e-11*RanGTP_cy + 1.2000000000000001e-11*RanGTP_nuc;
            break;
        case 1:
            dwdp[1] = -1.2000000000000001e-11*GDP*RCC1_Ran;
            break;
        case 2:
            dwdp[1] = 1.2000000000000001e-11*RCC1_RanGDP;
            break;
        case 3:
            dwdp[2] = -1.2000000000000001e-11*RCC1_RanGTP;
            break;
        case 4:
            dwdp[2] = 1.2000000000000001e-11*GTP*RCC1_Ran;
            break;
        case 5:
            dwdp[3] = -1.2000000000000001e-11*RanGDP_cy + 1.2000000000000001e-11*RanGDP_nuc;
            break;
        case 6:
            dwdp[4] = -1.2000000000000001e-11*RCC1_RanGDP;
            break;
        case 7:
            dwdp[4] = 1.2000000000000001e-11*RCC1*RanGDP_nuc;
            break;
        case 8:
            dwdp[5] = -1.7999999999999999e-11*RanBP1_RanGDP_kcat*RanGAP*RanGTP_RanBP1/pow(RanBP1_RanGDP_Km + RanGTP_RanBP1, 2);
            break;
        case 9:
            dwdp[5] = 1.7999999999999999e-11*RanGAP*RanGTP_RanBP1/(RanBP1_RanGDP_Km + RanGTP_RanBP1);
            break;
        case 10:
            dwdp[6] = -1.7999999999999999e-11*RanGAP*RanGAP_RanGDP_kcat_GAP*RanGTP_cy/pow(RanGAP_RanGDP_Km_GAP + RanGTP_cy, 2);
            break;
        case 11:
            dwdp[6] = 1.7999999999999999e-11*RanGAP*RanGTP_cy/(RanGAP_RanGDP_Km_GAP + RanGTP_cy);
            break;
        case 12:
            dwdp[7] = -1.7999999999999999e-11*RanGTP_RanBP1;
            break;
        case 13:
            dwdp[7] = 1.7999999999999999e-11*RanBP1*RanGTP_cy;
            break;
        case 14:
            dwdp[8] = -1.2000000000000001e-11*RCC1*RanGTP_nuc;
            break;
        case 15:
            dwdp[8] = 1.2000000000000001e-11*RCC1_RanGTP;
            break;
    }
}