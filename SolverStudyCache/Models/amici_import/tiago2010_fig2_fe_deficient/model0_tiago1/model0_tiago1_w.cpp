#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_tiago1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = re1_k1*s1;
    w[1] = re10_k1*s1;
    w[2] = re11_k1*s1;
    w[3] = re12_k1*s1;
    w[4] = re14_k1*s9;
    w[5] = re15_k1*s8;
    w[6] = re16_k1*s1;
    w[7] = re17_k1*s11;
    w[8] = re18_k1*s1;
    w[9] = re19_k1*s12;
    w[10] = re2_k1*s2;
    w[11] = re22_k1*s1;
    w[12] = re23_k1*s13;
    w[13] = re24_k1*s1;
    w[14] = re25_k1*s14;
    w[15] = re26_k1*s1;
    w[16] = re28_k1*s1;
    w[17] = re29_k1*s16;
    w[18] = re3_k1*s3;
    w[19] = re30_k1*s1;
    w[20] = re31_k1*s17;
    w[21] = re33_k1*s7;
    w[22] = re34_k1*s15;
    w[23] = re4_k1*s4;
    w[24] = re5_k1*s2;
    w[25] = re6_k1*s1;
    w[26] = re7_k1*s5;
    w[27] = re8_k1*s1;
    w[28] = re9_k1*s6;
}