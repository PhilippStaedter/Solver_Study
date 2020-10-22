#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_xie1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*PT*bccpt;
    dwdx[1] = 1.0*dcc;
    dwdx[2] = 1.0*ubcc;
    dwdx[3] = 1.0*bccperp*(-prcper + 1);
    dwdx[4] = 1.0*bccvrip*(-prcv + 1);
    dwdx[5] = 1.0*bccpdpp*(-prcpdp + 1);
    dwdx[6] = 1.0*bcctimp*(-prct + 1);
    dwdx[7] = 1.0*dccpt;
    dwdx[8] = 1.0*ubccpt;
    dwdx[9] = 1.0*CYC*bcc;
    dwdx[10] = 1.0*dclk;
    dwdx[11] = 1.0*CLK*bcc;
    dwdx[12] = 1.0*dpdp;
    dwdx[13] = 1.0*bpdpclkp*(-prpc - prvc + 1);
    dwdx[14] = 1.0*TIM*bpt;
    dwdx[15] = 1.0*dper;
    dwdx[16] = 1.0*CC*bccpt;
    dwdx[17] = 1.0*dpt;
    dwdx[18] = 1.0*ubpt;
    dwdx[19] = 1.0*PER*bpt;
    dwdx[20] = 1.0*dtim;
    dwdx[21] = 1.0*dvri;
    dwdx[22] = 1.0*bvriclkp*(-prpc - prvc + 1);
    dwdx[23] = 1.0*tlclk;
    dwdx[24] = 1.0*dclkm;
    dwdx[25] = 1.0*prpc*tcpdpclkp + 1.0*prvc*tcvriclkp + 1.0*tcclkp*(-prpc - prvc + 1);
    dwdx[26] = 1.0*dpdpm;
    dwdx[27] = 1.0*tlpdp;
    dwdx[28] = 1.0*tcccpdpp*(-pow(-prcpdp + 1, npdp) + 1) + 1.0*tcdvpmt*pow(-prcpdp + 1, npdp);
    dwdx[29] = 1.0*dperm;
    dwdx[30] = 1.0*tlper;
    dwdx[31] = 1.0*tcccperp*(-pow(-prcper + 1, npt) + 1) + 1.0*tcdvpmt*pow(-prcper + 1, npt);
    dwdx[32] = 1.0*pdpp*(npdp*tcccpdpp*pow(-prcpdp + 1, npdp)/(-prcpdp + 1) - npdp*tcdvpmt*pow(-prcpdp + 1, npdp)/(-prcpdp + 1));
    dwdx[33] = 1.0*ubccpdpp;
    dwdx[34] = -1.0*CC*bccpdpp;
    dwdx[35] = 1.0*perp*(npt*tcccperp*pow(-prcper + 1, npt)/(-prcper + 1) - npt*tcdvpmt*pow(-prcper + 1, npt)/(-prcper + 1));
    dwdx[36] = -1.0*CC*bccperp;
    dwdx[37] = 1.0*ubccperp;
    dwdx[38] = 1.0*timp*(npt*tccctimp*pow(-prct + 1, npt)/(-prct + 1) - npt*tcdvpmt*pow(-prct + 1, npt)/(-prct + 1));
    dwdx[39] = -1.0*CC*bcctimp;
    dwdx[40] = 1.0*ubcctimp;
    dwdx[41] = 1.0*vrip*(nvri*tcccvrip*pow(-prcv + 1, nvri)/(-prcv + 1) - nvri*tcdvpmt*pow(-prcv + 1, nvri)/(-prcv + 1));
    dwdx[42] = 1.0*ubccvrip;
    dwdx[43] = -1.0*CC*bccvrip;
    dwdx[44] = 1.0*clkp*(-tcclkp + tcpdpclkp);
    dwdx[45] = -1.0*VRI*bvriclkp;
    dwdx[46] = -1.0*PDP*bpdpclkp;
    dwdx[47] = 1.0*ubpdpclkp;
    dwdx[48] = 1.0*clkp*(-tcclkp + tcvriclkp);
    dwdx[49] = -1.0*VRI*bvriclkp;
    dwdx[50] = 1.0*ubvriclkp;
    dwdx[51] = -1.0*PDP*bpdpclkp;
    dwdx[52] = 1.0*dtimm;
    dwdx[53] = 1.0*tltim;
    dwdx[54] = 1.0*tccctimp*(-pow(-prct + 1, npt) + 1) + 1.0*tcdvpmt*pow(-prct + 1, npt);
    dwdx[55] = 1.0*dvrim;
    dwdx[56] = 1.0*tlvri;
    dwdx[57] = 1.0*tcccvrip*(-pow(-prcv + 1, nvri) + 1) + 1.0*tcdvpmt*pow(-prcv + 1, nvri);
}