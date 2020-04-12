#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_odea1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = IKK;
    x_solver[1] = IkBaIKK;
    x_solver[2] = IkBaIKKNFkB;
    x_solver[3] = IkBaNFkB_cytoplasm;
    x_solver[4] = IkBaNFkB_nucleus;
    x_solver[5] = IkBa_cytoplasm;
    x_solver[6] = IkBa_mRNA;
    x_solver[7] = IkBa_nucleus;
    x_solver[8] = IkBbIKK;
    x_solver[9] = IkBbIKKNFkB;
    x_solver[10] = IkBbNFkB_cytoplasm;
    x_solver[11] = IkBbNFkB_nucleus;
    x_solver[12] = IkBb_cytoplasm;
    x_solver[13] = IkBb_mRNA;
    x_solver[14] = IkBb_nucleus;
    x_solver[15] = IkBeIKK;
    x_solver[16] = IkBeIKKNFkB;
    x_solver[17] = IkBeNFkB_cytoplasm;
    x_solver[18] = IkBeNFkB_nucleus;
    x_solver[19] = IkBe_cytoplasm;
    x_solver[20] = IkBe_mRNA;
    x_solver[21] = IkBe_nucleus;
    x_solver[22] = NFkB_cytoplasm;
    x_solver[23] = NFkB_nucleus;
}