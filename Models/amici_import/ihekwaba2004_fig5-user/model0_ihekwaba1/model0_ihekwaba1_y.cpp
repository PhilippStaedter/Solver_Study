#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_ihekwaba1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = IKK;
    y[1] = IKKIkBa;
    y[2] = IKKIkBaNFkB;
    y[3] = IKKIkBb;
    y[4] = IKKIkBbNFkB;
    y[5] = IKKIkBe;
    y[6] = IKKIkBeNFkB;
    y[7] = IkBa;
    y[8] = IkBaNFkB;
    y[9] = IkBan;
    y[10] = IkBanNFkBn;
    y[11] = IkBat;
    y[12] = IkBb;
    y[13] = IkBbNFkB;
    y[14] = IkBbn;
    y[15] = IkBbnNFkBn;
    y[16] = IkBbt;
    y[17] = IkBe;
    y[18] = IkBeNFkB;
    y[19] = IkBen;
    y[20] = IkBenNFkBn;
    y[21] = IkBet;
    y[22] = NFkB;
    y[23] = NFkBn;
    y[24] = sink;
    y[25] = source;
}