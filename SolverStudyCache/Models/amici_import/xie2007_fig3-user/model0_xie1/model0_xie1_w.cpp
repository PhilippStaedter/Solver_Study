#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_xie1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*CC*PT*bccpt;
    w[1] = 1.0*clkm*tlclk;
    w[2] = 1.0*dperm*perm;
    w[3] = 1.0*perm*tlper;
    w[4] = 1.0*CC*dcc;
    w[5] = 1.0*perp*(tcccperp*(-pow(-prcper + 1, npt) + 1) + tcdvpmt*pow(-prcper + 1, npt));
    w[6] = 1.0*vrip*(tcccvrip*(-pow(-prcv + 1, nvri) + 1) + tcdvpmt*pow(-prcv + 1, nvri));
    w[7] = 1.0*CCPT*dccpt;
    w[8] = 1.0*dvrim*vrim;
    w[9] = 1.0*tlvri*vrim;
    w[10] = 1.0*VRI*dvri;
    w[11] = 1.0*pdpp*(tcccpdpp*(-pow(-prcpdp + 1, npdp) + 1) + tcdvpmt*pow(-prcpdp + 1, npdp));
    w[12] = 1.0*dpdpm*pdpm;
    w[13] = 1.0*pdpm*tlpdp;
    w[14] = 1.0*PDP*dpdp;
    w[15] = 1.0*PT*dpt;
    w[16] = 1.0*clkp*(prpc*tcpdpclkp + prvc*tcvriclkp + tcclkp*(-prpc - prvc + 1));
    w[17] = 1.0*CLK*CYC*bcc;
    w[18] = 1.0*CLK*dclk;
    w[19] = 1.0*CC*ubcc;
    w[20] = 1.0*PER*TIM*bpt;
    w[21] = 1.0*PT*ubpt;
    w[22] = 1.0*PER*dper;
    w[23] = 1.0*timp*(tccctimp*(-pow(-prct + 1, npt) + 1) + tcdvpmt*pow(-prct + 1, npt));
    w[24] = 1.0*dtimm*timm;
    w[25] = 1.0*timm*tltim;
    w[26] = 1.0*TIM*dtim;
    w[27] = 1.0*CCPT*ubccpt;
    w[28] = 1.0*CC*bccperp*(-prcper + 1);
    w[29] = 1.0*prcper*ubccperp;
    w[30] = 1.0*prcv*ubccvrip;
    w[31] = 1.0*CC*bccvrip*(-prcv + 1);
    w[32] = 1.0*prcpdp*ubccpdpp;
    w[33] = 1.0*CC*bccpdpp*(-prcpdp + 1);
    w[34] = 1.0*VRI*bvriclkp*(-prpc - prvc + 1);
    w[35] = 1.0*prvc*ubvriclkp;
    w[36] = 1.0*PDP*bpdpclkp*(-prpc - prvc + 1);
    w[37] = 1.0*prpc*ubpdpclkp;
    w[38] = 1.0*CC*bcctimp*(-prct + 1);
    w[39] = 1.0*prct*ubcctimp;
    w[40] = 1.0*clkm*dclkm;
}