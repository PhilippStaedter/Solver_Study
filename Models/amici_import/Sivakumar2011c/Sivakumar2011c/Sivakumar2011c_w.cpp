#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Sivakumar2011c(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = kass_r68*s305;
    w[1] = kass_r1*s1*s5 - kdiss_r1*s16;
    w[2] = kass_r5*s16*s28 - kdiss_r5*s27;
    w[3] = kass_r47*s121*s36 - kdiss_r47*s123;
    w[4] = kass_r48*s123*s46 - kdiss_r48*s129;
    w[5] = kass_r54*s123*s75 - kdiss_r54*s159;
    w[6] = kass_r58*s36 - kdiss_r58*s232;
    w[7] = kass_r63*s174*s232 - kdiss_r63*s176;
    w[8] = kass_r66*s173*s183 - kdiss_r66*s188;
    w[9] = kass_r88*s252*s61 - kdiss_r88*s259;
    w[10] = kass_r90*s259*s268 - kdiss_r90*s266;
    w[11] = kass_r91*s266 - kdiss_r91*s155*s267;
    w[12] = kass_r92*s267 - kdiss_r92*s260*s61;
    w[13] = kass_r96*s159*s268 - kdiss_r96*s275;
    w[14] = kass_r98*s275 - kdiss_r98*s101*s278;
    w[15] = kass_r99*s278 - kdiss_r99*s164*s270;
    w[16] = kass_r102*s286*s31 - kdiss_r102*s288;
    w[17] = kass_r103*s102*s288 - kdiss_r103*s292;
    w[18] = kass_r105*s292 - kdiss_r105*s37;
    w[19] = kass_r106*s286 - kdiss_r106*s30;
    w[20] = kass_r107*s239 - kdiss_r107*s5;
    w[21] = s30*(kass_r104_s30*s107*s32 - kdiss_r104_s30*s286*s33);
    w[22] = s30*(kass_r85_s30*s129*s32 - kdiss_r85_s30*s245*s33);
    w[23] = kass_r65*s171*s179 - kdiss_r65*s183;
    w[24] = kass_r64*s170*s176 - kdiss_r64*s179;
    w[25] = kass_re65*s260;
    w[26] = kass_re64*s270;
    w[27] = kass_r67*s172*s188 - kdiss_r67*s305;
    w[28] = kI_r86_s304*s37*(kass_r86_s37*s245*pow(s32, 3) - kdiss_r86_s37*s252*pow(s33, 3))/(kI_r86_s304 + s304);
}