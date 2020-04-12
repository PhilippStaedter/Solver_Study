#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_ihekwaba1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = IKK;
    x_solver[1] = IKKIkBa;
    x_solver[2] = IKKIkBaNFkB;
    x_solver[3] = IKKIkBb;
    x_solver[4] = IKKIkBbNFkB;
    x_solver[5] = IKKIkBe;
    x_solver[6] = IKKIkBeNFkB;
    x_solver[7] = IkBa;
    x_solver[8] = IkBaNFkB;
    x_solver[9] = IkBan;
    x_solver[10] = IkBanNFkBn;
    x_solver[11] = IkBat;
    x_solver[12] = IkBb;
    x_solver[13] = IkBbNFkB;
    x_solver[14] = IkBbn;
    x_solver[15] = IkBbnNFkBn;
    x_solver[16] = IkBbt;
    x_solver[17] = IkBe;
    x_solver[18] = IkBeNFkB;
    x_solver[19] = IkBen;
    x_solver[20] = IkBenNFkBn;
    x_solver[21] = IkBet;
    x_solver[22] = NFkB;
    x_solver[23] = NFkBn;
    x_solver[24] = sink;
    x_solver[25] = source;
}