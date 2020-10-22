#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Holzhutter2004(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Glcin;
    x_rdata[1] = MgATP;
    x_rdata[2] = Glc6P;
    x_rdata[3] = MgADP;
    x_rdata[4] = Fru6P;
    x_rdata[5] = Fru16P2;
    x_rdata[6] = GraP;
    x_rdata[7] = DHAP;
    x_rdata[8] = Phi;
    x_rdata[9] = NAD;
    x_rdata[10] = Gri13P2;
    x_rdata[11] = NADH;
    x_rdata[12] = Gri3P;
    x_rdata[13] = Gri23P2f;
    x_rdata[14] = Gri2P;
    x_rdata[15] = PEP;
    x_rdata[16] = Pyr;
    x_rdata[17] = Lac;
    x_rdata[18] = NADPHf;
    x_rdata[19] = NADPf;
    x_rdata[20] = AMPf;
    x_rdata[21] = ADPf;
    x_rdata[22] = GlcA6P;
    x_rdata[23] = Rul5P;
    x_rdata[24] = GSSG;
    x_rdata[25] = GSH;
    x_rdata[26] = Xul5P;
    x_rdata[27] = Rib5P;
    x_rdata[28] = Sed7P;
    x_rdata[29] = E4P;
    x_rdata[30] = MgAMP;
    x_rdata[31] = ATPf;
    x_rdata[32] = Mgf;
    x_rdata[33] = MgGri23P2;
    x_rdata[34] = P1NADP;
    x_rdata[35] = P1f;
    x_rdata[36] = P1NADPH;
    x_rdata[37] = P2NADP;
    x_rdata[38] = P2f;
    x_rdata[39] = P2NADPH;
    x_rdata[40] = PRPP;
    x_rdata[41] = Lacex;
    x_rdata[42] = Pyrex;
    x_rdata[43] = Glcout;
    x_rdata[44] = Phiex;
}