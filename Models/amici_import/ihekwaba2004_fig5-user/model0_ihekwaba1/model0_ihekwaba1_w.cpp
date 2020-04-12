#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_ihekwaba1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = IKKIkBa*NFkB*k7;
    w[1] = IKK*IkBb*k40;
    w[2] = IKKIkBb*k63;
    w[3] = IKKIkBb*k41;
    w[4] = IKKIkBe*NFkB*k13;
    w[5] = IKKIkBeNFkB*k15;
    w[6] = IKKIkBeNFkB*k14;
    w[7] = IKK*IkBe*k46;
    w[8] = IKKIkBe*k64;
    w[9] = IKKIkBe*k47;
    w[10] = IKK*IkBaNFkB*k52;
    w[11] = IKKIkBaNFkB*k9;
    w[12] = IKKIkBaNFkB*k53;
    w[13] = IKK*IkBbNFkB*k55;
    w[14] = IKKIkBbNFkB*k56;
    w[15] = IKK*IkBeNFkB*k58;
    w[16] = IKKIkBeNFkB*k59;
    w[17] = IkBa*NFkB*k1;
    w[18] = pow(NFkBn, 2)*k28;
    w[19] = IkBaNFkB*k2;
    w[20] = IkBan*NFkBn*k21;
    w[21] = IkBanNFkBn*k22;
    w[22] = IKKIkBaNFkB*k8;
    w[23] = IkBanNFkBn*k54;
    w[24] = IkBaNFkB*k16;
    w[25] = IkBat*k29;
    w[26] = IkBan*k39;
    w[27] = IkBa*k38;
    w[28] = IkBat*k36;
    w[29] = IkBb*NFkB*k3;
    w[30] = IkBbNFkB*k4;
    w[31] = IkBbn*NFkBn*k23;
    w[32] = IkBbnNFkBn*k24;
    w[33] = IKK*IkBa*k34;
    w[34] = IkBbnNFkBn*k57;
    w[35] = IkBbNFkB*k17;
    w[36] = IkBbt*k31;
    w[37] = IkBbn*k45;
    w[38] = IkBb*k44;
    w[39] = IkBbt*k42;
    w[40] = IkBe*NFkB*k5;
    w[41] = IkBeNFkB*k6;
    w[42] = IkBen*NFkBn*k25;
    w[43] = IkBenNFkBn*k26;
    w[44] = IKKIkBa*k62;
    w[45] = IkBenNFkBn*k60;
    w[46] = IkBeNFkB*k18;
    w[47] = IkBet*k33;
    w[48] = IkBen*k51;
    w[49] = IkBe*k50;
    w[50] = IkBet*k48;
    w[51] = NFkBn*k20;
    w[52] = NFkB*k19;
    w[53] = k27*source;
    w[54] = k30*source;
    w[55] = IKKIkBa*k35;
    w[56] = k32*source;
    w[57] = IkBa*k37;
    w[58] = IkBb*k43;
    w[59] = IkBe*k49;
    w[60] = IKK*k61;
    w[61] = IKKIkBb*NFkB*k10;
    w[62] = IKKIkBbNFkB*k12;
    w[63] = IKKIkBbNFkB*k11;
}