#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_tiago2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = s1;
            break;
        case 1:
            dwdp[1] = s1;
            break;
        case 2:
            dwdp[2] = s1;
            break;
        case 3:
            dwdp[3] = s1;
            break;
        case 4:
            dwdp[4] = s9;
            break;
        case 5:
            dwdp[5] = s8;
            break;
        case 6:
            dwdp[6] = s1;
            break;
        case 7:
            dwdp[7] = s11;
            break;
        case 8:
            dwdp[8] = s1;
            break;
        case 9:
            dwdp[9] = s12;
            break;
        case 10:
            dwdp[10] = s2;
            break;
        case 11:
            dwdp[11] = s1;
            break;
        case 12:
            dwdp[12] = s13;
            break;
        case 13:
            dwdp[13] = s1;
            break;
        case 14:
            dwdp[14] = s14;
            break;
        case 15:
            dwdp[15] = s1;
            break;
        case 16:
            dwdp[16] = s1;
            break;
        case 17:
            dwdp[17] = s16;
            break;
        case 18:
            dwdp[18] = s3;
            break;
        case 19:
            dwdp[19] = s1;
            break;
        case 20:
            dwdp[20] = s17;
            break;
        case 21:
            dwdp[21] = s7;
            break;
        case 22:
            dwdp[22] = s15;
            break;
        case 23:
            dwdp[23] = s4;
            break;
        case 24:
            dwdp[24] = s2;
            break;
        case 25:
            dwdp[25] = s1;
            break;
        case 26:
            dwdp[26] = s5;
            break;
        case 27:
            dwdp[27] = s1;
            break;
        case 28:
            dwdp[28] = s6;
            break;
    }
}