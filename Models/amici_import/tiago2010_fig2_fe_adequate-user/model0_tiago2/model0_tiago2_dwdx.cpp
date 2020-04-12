#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_tiago2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = re1_k1;
    dwdx[1] = re10_k1;
    dwdx[2] = re11_k1;
    dwdx[3] = re12_k1;
    dwdx[4] = re16_k1;
    dwdx[5] = re18_k1;
    dwdx[6] = re22_k1;
    dwdx[7] = re24_k1;
    dwdx[8] = re26_k1;
    dwdx[9] = re28_k1;
    dwdx[10] = re30_k1;
    dwdx[11] = re6_k1;
    dwdx[12] = re8_k1;
    dwdx[13] = re17_k1;
    dwdx[14] = re19_k1;
    dwdx[15] = re23_k1;
    dwdx[16] = re25_k1;
    dwdx[17] = re33_k1;
    dwdx[18] = re29_k1;
    dwdx[19] = re31_k1;
    dwdx[20] = re2_k1;
    dwdx[21] = re5_k1;
    dwdx[22] = re3_k1;
    dwdx[23] = re4_k1;
    dwdx[24] = re7_k1;
    dwdx[25] = re9_k1;
    dwdx[26] = re32_k1;
    dwdx[27] = re15_k1;
    dwdx[28] = re14_k1;
}