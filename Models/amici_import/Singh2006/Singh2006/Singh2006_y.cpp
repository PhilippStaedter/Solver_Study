#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Singh2006(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = x1;
    y[1] = x2;
    y[2] = x6;
    y[3] = x5;
    y[4] = x7;
    y[5] = x8;
    y[6] = x16;
    y[7] = x15;
    y[8] = x9;
    y[9] = x11;
    y[10] = x10;
    y[11] = x12;
    y[12] = x29;
    y[13] = x30;
    y[14] = x39;
    y[15] = x46;
    y[16] = x40;
    y[17] = x45;
    y[18] = x41;
    y[19] = x44;
    y[20] = x18;
    y[21] = x17;
    y[22] = x14;
    y[23] = x22;
    y[24] = x32;
    y[25] = x13;
    y[26] = x20;
    y[27] = x21;
    y[28] = x23;
    y[29] = x27;
    y[30] = x24;
    y[31] = x25;
    y[32] = x34;
    y[33] = x36;
    y[34] = x42;
    y[35] = x37;
    y[36] = x47;
    y[37] = x48;
    y[38] = x52;
    y[39] = x51;
    y[40] = x50;
    y[41] = x53;
    y[42] = x54;
    y[43] = x55;
    y[44] = x60;
    y[45] = x59;
    y[46] = x57;
    y[47] = x61;
    y[48] = x62;
    y[49] = x63;
    y[50] = x65;
    y[51] = x68;
    y[52] = x35;
    y[53] = x28;
    y[54] = x31;
    y[55] = x56;
    y[56] = x43;
    y[57] = x3;
    y[58] = x58;
    y[59] = x4;
    y[60] = x26;
    y[61] = x49;
    y[62] = x33;
    y[63] = x64;
    y[64] = x19;
    y[65] = x38;
    y[66] = x66;
    y[67] = x67;
}