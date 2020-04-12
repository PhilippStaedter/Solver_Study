#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_odea1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = IKK;
    x_rdata[1] = IkBaIKK;
    x_rdata[2] = IkBaIKKNFkB;
    x_rdata[3] = IkBaNFkB_cytoplasm;
    x_rdata[4] = IkBaNFkB_nucleus;
    x_rdata[5] = IkBa_cytoplasm;
    x_rdata[6] = IkBa_mRNA;
    x_rdata[7] = IkBa_nucleus;
    x_rdata[8] = IkBbIKK;
    x_rdata[9] = IkBbIKKNFkB;
    x_rdata[10] = IkBbNFkB_cytoplasm;
    x_rdata[11] = IkBbNFkB_nucleus;
    x_rdata[12] = IkBb_cytoplasm;
    x_rdata[13] = IkBb_mRNA;
    x_rdata[14] = IkBb_nucleus;
    x_rdata[15] = IkBeIKK;
    x_rdata[16] = IkBeIKKNFkB;
    x_rdata[17] = IkBeNFkB_cytoplasm;
    x_rdata[18] = IkBeNFkB_nucleus;
    x_rdata[19] = IkBe_cytoplasm;
    x_rdata[20] = IkBe_mRNA;
    x_rdata[21] = IkBe_nucleus;
    x_rdata[22] = NFkB_cytoplasm;
    x_rdata[23] = NFkB_nucleus;
}