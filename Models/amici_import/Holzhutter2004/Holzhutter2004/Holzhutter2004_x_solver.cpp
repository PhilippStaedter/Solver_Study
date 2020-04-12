#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Holzhutter2004(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Glcin;
    x_solver[1] = MgATP;
    x_solver[2] = Glc6P;
    x_solver[3] = MgADP;
    x_solver[4] = Fru6P;
    x_solver[5] = Fru16P2;
    x_solver[6] = GraP;
    x_solver[7] = DHAP;
    x_solver[8] = Phi;
    x_solver[9] = NAD;
    x_solver[10] = Gri13P2;
    x_solver[11] = NADH;
    x_solver[12] = Gri3P;
    x_solver[13] = Gri23P2f;
    x_solver[14] = Gri2P;
    x_solver[15] = PEP;
    x_solver[16] = Pyr;
    x_solver[17] = Lac;
    x_solver[18] = NADPHf;
    x_solver[19] = NADPf;
    x_solver[20] = AMPf;
    x_solver[21] = ADPf;
    x_solver[22] = GlcA6P;
    x_solver[23] = Rul5P;
    x_solver[24] = GSSG;
    x_solver[25] = GSH;
    x_solver[26] = Xul5P;
    x_solver[27] = Rib5P;
    x_solver[28] = Sed7P;
    x_solver[29] = E4P;
    x_solver[30] = MgAMP;
    x_solver[31] = ATPf;
    x_solver[32] = Mgf;
    x_solver[33] = MgGri23P2;
    x_solver[34] = P1NADP;
    x_solver[35] = P1f;
    x_solver[36] = P1NADPH;
    x_solver[37] = P2NADP;
    x_solver[38] = P2f;
    x_solver[39] = P2NADPH;
    x_solver[40] = PRPP;
    x_solver[41] = Lacex;
    x_solver[42] = Pyrex;
    x_solver[43] = Glcout;
    x_solver[44] = Phiex;
}