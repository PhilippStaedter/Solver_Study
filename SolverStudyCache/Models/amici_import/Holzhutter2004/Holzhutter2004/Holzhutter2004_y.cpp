#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Holzhutter2004(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Glcin;
    y[1] = MgATP;
    y[2] = Glc6P;
    y[3] = MgADP;
    y[4] = Fru6P;
    y[5] = Fru16P2;
    y[6] = GraP;
    y[7] = DHAP;
    y[8] = Phi;
    y[9] = NAD;
    y[10] = Gri13P2;
    y[11] = NADH;
    y[12] = Gri3P;
    y[13] = Gri23P2f;
    y[14] = Gri2P;
    y[15] = PEP;
    y[16] = Pyr;
    y[17] = Lac;
    y[18] = NADPHf;
    y[19] = NADPf;
    y[20] = AMPf;
    y[21] = ADPf;
    y[22] = GlcA6P;
    y[23] = Rul5P;
    y[24] = GSSG;
    y[25] = GSH;
    y[26] = Xul5P;
    y[27] = Rib5P;
    y[28] = Sed7P;
    y[29] = E4P;
    y[30] = MgAMP;
    y[31] = ATPf;
    y[32] = Mgf;
    y[33] = MgGri23P2;
    y[34] = P1NADP;
    y[35] = P1f;
    y[36] = P1NADPH;
    y[37] = P2NADP;
    y[38] = P2f;
    y[39] = P2NADPH;
    y[40] = PRPP;
    y[41] = Lacex;
    y[42] = Pyrex;
    y[43] = Glcout;
    y[44] = Phiex;
}