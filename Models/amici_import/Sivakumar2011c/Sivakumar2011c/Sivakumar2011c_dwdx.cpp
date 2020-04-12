#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Sivakumar2011c(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = kass_r1*s1;
    dwdx[1] = -kdiss_r107;
    dwdx[2] = -kdiss_r1;
    dwdx[3] = kass_r5*s28;
    dwdx[4] = kass_r1*s5;
    dwdx[5] = -kdiss_r5;
    dwdx[6] = kass_r5*s16;
    dwdx[7] = -kdiss_r106;
    dwdx[8] = kass_r104_s30*s107*s32 - kdiss_r104_s30*s286*s33;
    dwdx[9] = kass_r85_s30*s129*s32 - kdiss_r85_s30*s245*s33;
    dwdx[10] = kass_r102*s286;
    dwdx[11] = kass_r104_s30*s107*s30;
    dwdx[12] = kass_r85_s30*s129*s30;
    dwdx[13] = 3*kI_r86_s304*kass_r86_s37*s245*pow(s32, 2)*s37/(kI_r86_s304 + s304);
    dwdx[14] = -kdiss_r104_s30*s286*s30;
    dwdx[15] = -kdiss_r85_s30*s245*s30;
    dwdx[16] = -3*kI_r86_s304*kdiss_r86_s37*s252*pow(s33, 2)*s37/(kI_r86_s304 + s304);
    dwdx[17] = -kdiss_r105;
    dwdx[18] = kI_r86_s304*(kass_r86_s37*s245*pow(s32, 3) - kdiss_r86_s37*s252*pow(s33, 3))/(kI_r86_s304 + s304);
    dwdx[19] = kass_r48*s123;
    dwdx[20] = kass_r54*s123;
    dwdx[21] = -kdiss_r98*s278;
    dwdx[22] = kass_r103*s288;
    dwdx[23] = kass_r104_s30*s30*s32;
    dwdx[24] = kass_r47*s36;
    dwdx[25] = -kdiss_r91*s267;
    dwdx[26] = -kdiss_r99*s270;
    dwdx[27] = kass_r65*s179;
    dwdx[28] = kass_r67*s188;
    dwdx[29] = kass_r66*s183;
    dwdx[30] = kass_r64*s176;
    dwdx[31] = kass_r63*s232;
    dwdx[32] = kass_r107;
    dwdx[33] = kass_r47*s121;
    dwdx[34] = kass_r58;
    dwdx[35] = -kdiss_r47;
    dwdx[36] = kass_r48*s46;
    dwdx[37] = kass_r54*s75;
    dwdx[38] = -kdiss_r48;
    dwdx[39] = kass_r85_s30*s30*s32;
    dwdx[40] = -kdiss_r54;
    dwdx[41] = kass_r96*s268;
    dwdx[42] = -kdiss_r58;
    dwdx[43] = kass_r63*s174;
    dwdx[44] = -kdiss_r63;
    dwdx[45] = kass_r64*s170;
    dwdx[46] = kass_r65*s171;
    dwdx[47] = -kdiss_r64;
    dwdx[48] = kass_r66*s173;
    dwdx[49] = -kdiss_r65;
    dwdx[50] = -kdiss_r66;
    dwdx[51] = kass_r67*s172;
    dwdx[52] = -kdiss_r85_s30*s30*s33;
    dwdx[53] = kI_r86_s304*kass_r86_s37*pow(s32, 3)*s37/(kI_r86_s304 + s304);
    dwdx[54] = kass_r88*s61;
    dwdx[55] = -kI_r86_s304*kdiss_r86_s37*pow(s33, 3)*s37/(kI_r86_s304 + s304);
    dwdx[56] = kass_r90*s259;
    dwdx[57] = kass_r96*s159;
    dwdx[58] = -kdiss_r92*s61;
    dwdx[59] = kass_re65;
    dwdx[60] = -kdiss_r99*s164;
    dwdx[61] = kass_re64;
    dwdx[62] = -kdiss_r96;
    dwdx[63] = kass_r98;
    dwdx[64] = -kdiss_r98*s101;
    dwdx[65] = kass_r99;
    dwdx[66] = kass_r102*s31;
    dwdx[67] = kass_r106;
    dwdx[68] = -kdiss_r104_s30*s30*s33;
    dwdx[69] = -kdiss_r102;
    dwdx[70] = kass_r103*s102;
    dwdx[71] = -kdiss_r103;
    dwdx[72] = kass_r105;
    dwdx[73] = kass_r88*s252;
    dwdx[74] = -kdiss_r92*s260;
    dwdx[75] = -kdiss_r88;
    dwdx[76] = kass_r90*s268;
    dwdx[77] = -kdiss_r90;
    dwdx[78] = kass_r91;
    dwdx[79] = -kdiss_r91*s155;
    dwdx[80] = kass_r92;
    dwdx[81] = -kI_r86_s304*s37*(kass_r86_s37*s245*pow(s32, 3) - kdiss_r86_s37*s252*pow(s33, 3))/pow(kI_r86_s304 + s304, 2);
    dwdx[82] = kass_r68;
    dwdx[83] = -kdiss_r67;
}