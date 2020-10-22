#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_ihekwaba1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = IKK;
    x_rdata[1] = IKKIkBa;
    x_rdata[2] = IKKIkBaNFkB;
    x_rdata[3] = IKKIkBb;
    x_rdata[4] = IKKIkBbNFkB;
    x_rdata[5] = IKKIkBe;
    x_rdata[6] = IKKIkBeNFkB;
    x_rdata[7] = IkBa;
    x_rdata[8] = IkBaNFkB;
    x_rdata[9] = IkBan;
    x_rdata[10] = IkBanNFkBn;
    x_rdata[11] = IkBat;
    x_rdata[12] = IkBb;
    x_rdata[13] = IkBbNFkB;
    x_rdata[14] = IkBbn;
    x_rdata[15] = IkBbnNFkBn;
    x_rdata[16] = IkBbt;
    x_rdata[17] = IkBe;
    x_rdata[18] = IkBeNFkB;
    x_rdata[19] = IkBen;
    x_rdata[20] = IkBenNFkBn;
    x_rdata[21] = IkBet;
    x_rdata[22] = NFkB;
    x_rdata[23] = NFkBn;
    x_rdata[24] = sink;
    x_rdata[25] = source;
}