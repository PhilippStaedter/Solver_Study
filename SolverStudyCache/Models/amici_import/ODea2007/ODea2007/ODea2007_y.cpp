#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_ODea2007(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = IkBa_mRNA;
    y[1] = IkBa_cytoplasm;
    y[2] = IkBa_nucleus;
    y[3] = IkBaIKK;
    y[4] = IkBaNFkB_cytoplasm;
    y[5] = IkBaNFkB_nucleus;
    y[6] = IkBaIKKNFkB;
    y[7] = NFkB_cytoplasm;
    y[8] = IKK;
    y[9] = NFkB_nucleus;
    y[10] = IkBbIKK;
    y[11] = IkBbIKKNFkB;
    y[12] = IkBbNFkB_nucleus;
    y[13] = IkBbNFkB_cytoplasm;
    y[14] = IkBb_nucleus;
    y[15] = IkBb_cytoplasm;
    y[16] = IkBb_mRNA;
    y[17] = IkBe_mRNA;
    y[18] = IkBe_cytoplasm;
    y[19] = IkBe_nucleus;
    y[20] = IkBeNFkB_cytoplasm;
    y[21] = IkBeNFkB_nucleus;
    y[22] = IkBeIKKNFkB;
    y[23] = IkBeIKK;
}