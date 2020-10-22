#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_ihekwaba1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = IkBb*k40;
    dwdx[1] = IkBe*k46;
    dwdx[2] = IkBaNFkB*k52;
    dwdx[3] = IkBbNFkB*k55;
    dwdx[4] = IkBeNFkB*k58;
    dwdx[5] = IkBa*k34;
    dwdx[6] = k61;
    dwdx[7] = NFkB*k7;
    dwdx[8] = k62;
    dwdx[9] = k35;
    dwdx[10] = k9;
    dwdx[11] = k53;
    dwdx[12] = k8;
    dwdx[13] = k63;
    dwdx[14] = k41;
    dwdx[15] = NFkB*k10;
    dwdx[16] = k56;
    dwdx[17] = k12;
    dwdx[18] = k11;
    dwdx[19] = NFkB*k13;
    dwdx[20] = k64;
    dwdx[21] = k47;
    dwdx[22] = k15;
    dwdx[23] = k14;
    dwdx[24] = k59;
    dwdx[25] = NFkB*k1;
    dwdx[26] = k38;
    dwdx[27] = IKK*k34;
    dwdx[28] = k37;
    dwdx[29] = IKK*k52;
    dwdx[30] = k2;
    dwdx[31] = k16;
    dwdx[32] = NFkBn*k21;
    dwdx[33] = k39;
    dwdx[34] = k22;
    dwdx[35] = k54;
    dwdx[36] = k29;
    dwdx[37] = k36;
    dwdx[38] = IKK*k40;
    dwdx[39] = NFkB*k3;
    dwdx[40] = k44;
    dwdx[41] = k43;
    dwdx[42] = IKK*k55;
    dwdx[43] = k4;
    dwdx[44] = k17;
    dwdx[45] = NFkBn*k23;
    dwdx[46] = k45;
    dwdx[47] = k24;
    dwdx[48] = k57;
    dwdx[49] = k31;
    dwdx[50] = k42;
    dwdx[51] = IKK*k46;
    dwdx[52] = NFkB*k5;
    dwdx[53] = k50;
    dwdx[54] = k49;
    dwdx[55] = IKK*k58;
    dwdx[56] = k6;
    dwdx[57] = k18;
    dwdx[58] = NFkBn*k25;
    dwdx[59] = k51;
    dwdx[60] = k26;
    dwdx[61] = k60;
    dwdx[62] = k33;
    dwdx[63] = k48;
    dwdx[64] = IKKIkBa*k7;
    dwdx[65] = IKKIkBe*k13;
    dwdx[66] = IkBa*k1;
    dwdx[67] = IkBb*k3;
    dwdx[68] = IkBe*k5;
    dwdx[69] = k19;
    dwdx[70] = IKKIkBb*k10;
    dwdx[71] = 2*NFkBn*k28;
    dwdx[72] = IkBan*k21;
    dwdx[73] = IkBbn*k23;
    dwdx[74] = IkBen*k25;
    dwdx[75] = k20;
    dwdx[76] = k27;
    dwdx[77] = k30;
    dwdx[78] = k32;
}