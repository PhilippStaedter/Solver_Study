#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_odea1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = IKK;
    y[1] = IkBaIKK;
    y[2] = IkBaIKKNFkB;
    y[3] = IkBaNFkB_cytoplasm;
    y[4] = IkBaNFkB_nucleus;
    y[5] = IkBa_cytoplasm;
    y[6] = IkBa_mRNA;
    y[7] = IkBa_nucleus;
    y[8] = IkBbIKK;
    y[9] = IkBbIKKNFkB;
    y[10] = IkBbNFkB_cytoplasm;
    y[11] = IkBbNFkB_nucleus;
    y[12] = IkBb_cytoplasm;
    y[13] = IkBb_mRNA;
    y[14] = IkBb_nucleus;
    y[15] = IkBeIKK;
    y[16] = IkBeIKKNFkB;
    y[17] = IkBeNFkB_cytoplasm;
    y[18] = IkBeNFkB_nucleus;
    y[19] = IkBe_cytoplasm;
    y[20] = IkBe_mRNA;
    y[21] = IkBe_nucleus;
    y[22] = NFkB_cytoplasm;
    y[23] = NFkB_nucleus;
}