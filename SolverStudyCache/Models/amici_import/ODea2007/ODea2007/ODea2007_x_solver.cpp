#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_ODea2007(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = IkBa_mRNA;
    x_solver[1] = IkBa_cytoplasm;
    x_solver[2] = IkBa_nucleus;
    x_solver[3] = IkBaIKK;
    x_solver[4] = IkBaNFkB_cytoplasm;
    x_solver[5] = IkBaNFkB_nucleus;
    x_solver[6] = IkBaIKKNFkB;
    x_solver[7] = NFkB_cytoplasm;
    x_solver[8] = IKK;
    x_solver[9] = NFkB_nucleus;
    x_solver[10] = IkBbIKK;
    x_solver[11] = IkBbIKKNFkB;
    x_solver[12] = IkBbNFkB_nucleus;
    x_solver[13] = IkBbNFkB_cytoplasm;
    x_solver[14] = IkBb_nucleus;
    x_solver[15] = IkBb_cytoplasm;
    x_solver[16] = IkBb_mRNA;
    x_solver[17] = IkBe_mRNA;
    x_solver[18] = IkBe_cytoplasm;
    x_solver[19] = IkBe_nucleus;
    x_solver[20] = IkBeNFkB_cytoplasm;
    x_solver[21] = IkBeNFkB_nucleus;
    x_solver[22] = IkBeIKKNFkB;
    x_solver[23] = IkBeIKK;
}