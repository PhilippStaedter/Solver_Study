#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_ODea2007(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = IkBa_mRNA;
    x_rdata[1] = IkBa_cytoplasm;
    x_rdata[2] = IkBa_nucleus;
    x_rdata[3] = IkBaIKK;
    x_rdata[4] = IkBaNFkB_cytoplasm;
    x_rdata[5] = IkBaNFkB_nucleus;
    x_rdata[6] = IkBaIKKNFkB;
    x_rdata[7] = NFkB_cytoplasm;
    x_rdata[8] = IKK;
    x_rdata[9] = NFkB_nucleus;
    x_rdata[10] = IkBbIKK;
    x_rdata[11] = IkBbIKKNFkB;
    x_rdata[12] = IkBbNFkB_nucleus;
    x_rdata[13] = IkBbNFkB_cytoplasm;
    x_rdata[14] = IkBb_nucleus;
    x_rdata[15] = IkBb_cytoplasm;
    x_rdata[16] = IkBb_mRNA;
    x_rdata[17] = IkBe_mRNA;
    x_rdata[18] = IkBe_cytoplasm;
    x_rdata[19] = IkBe_nucleus;
    x_rdata[20] = IkBeNFkB_cytoplasm;
    x_rdata[21] = IkBeNFkB_nucleus;
    x_rdata[22] = IkBeIKKNFkB;
    x_rdata[23] = IkBeIKK;
}