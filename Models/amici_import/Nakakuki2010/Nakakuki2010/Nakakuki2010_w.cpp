#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Nakakuki2010(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 940.0*ERK_c*Fct*V1*pMEK/(ERK_c + K1*(1 + pERK_c/K2));
    w[1] = 940.0*Fct*V2*pERK_c*pMEK/(K2*(ERK_c/K1 + 1) + pERK_c);
    w[2] = 940.0*V3*pERK_c/(K3*(1 + ppERK_c/K4) + pERK_c);
    w[3] = 940.0*V4*ppERK_c/(K4*(1 + pERK_c/K3) + ppERK_c);
    w[4] = 220.0*V5*pERK_n/(K5*(1 + ppERK_n/K6) + pERK_n);
    w[5] = 220.0*V6*ppERK_n/(K6*(1 + pERK_n/K5) + ppERK_n);
    w[6] = ERK_c*KimERK*Vc - ERK_n*KexERK*Vn;
    w[7] = -KexERKP*Vn*pERK_n + KimERKP*Vc*pERK_c;
    w[8] = -KexERKPP*Vn*ppERK_n + KimERKPP*Vc*ppERK_c;
    w[9] = 220.0*V10*pow(ppERK_n, n10)/(pow(K10, n10) + pow(ppERK_n, n10));
    w[10] = PreDUSPmRNA*Vn*p11;
    w[11] = 940.0*DUSPmRNA*p12;
    w[12] = 940.0*DUSPmRNA*p13;
    w[13] = 940.0*DUSP_c*V14*ppERK_c/(DUSP_c + K14);
    w[14] = 940.0*V15*pDUSP_c/(K15 + pDUSP_c);
    w[15] = 940.0*DUSP_c*p16;
    w[16] = 940.0*p17*pDUSP_c;
    w[17] = DUSP_c*KimDUSP*Vc - DUSP_n*KexDUSP*Vn;
    w[18] = -KexDUSPP*Vn*pDUSP_n + KimDUSPP*Vc*pDUSP_c;
    w[19] = 220.0*DUSP_n*V20*ppERK_n/(DUSP_n + K20);
    w[20] = 220.0*V21*pDUSP_n/(K21 + pDUSP_n);
    w[21] = 220.0*DUSP_n*p22;
    w[22] = 220.0*p23*pDUSP_n;
    w[23] = 940.0*RSK_c*V24*ppERK_c/(K24 + RSK_c);
    w[24] = 940.0*V25*pRSK_c/(K25 + pRSK_c);
    w[25] = -KexRSKP*Vn*pRSK_n + KimRSKP*Vc*pRSK_c;
    w[26] = 220.0*CREB_n*V27*pRSK_n/(CREB_n + K27);
    w[27] = 220.0*V28*pCREB_n/(K28 + pCREB_n);
    w[28] = 220.0*Elk1_n*V29*ppERK_n/(Elk1_n + K29);
    w[29] = 220.0*V30*pElk1_n/(K30 + pElk1_n);
    w[30] = 220.0*V31*pow(pCREB_n*pElk1_n, n31)/(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31));
    w[31] = PreFOSmRNA*Vn*p32;
    w[32] = 940.0*c_FOSmRNA*p33;
    w[33] = 940.0*c_FOSmRNA*p34;
    w[34] = 940.0*V35*c_FOS_c*ppERK_c/(K35 + c_FOS_c);
    w[35] = 940.0*V36*c_FOS_c*pRSK_c/(K36 + c_FOS_c);
    w[36] = 940.0*V37*pc_FOS_c/(K37 + pc_FOS_c);
    w[37] = 940.0*c_FOS_c*p38;
    w[38] = 940.0*p39*pc_FOS_c;
    w[39] = -FOSn*KexFOS*Vn + KimFOS*Vc*c_FOS_c;
    w[40] = -FOSn_2*KexFOSP*Vn + KimFOSP*Vc*pc_FOS_c;
    w[41] = 220.0*FOSn*V42*ppERK_n/(FOSn + K42);
    w[42] = 220.0*FOSn*V43*pRSK_n/(FOSn + K43);
    w[43] = 220.0*FOSn_2*V44/(FOSn_2 + K44);
    w[44] = 220.0*FOSn*p45;
    w[45] = 220.0*FOSn_2*p46;
    w[46] = -220.0*m52*pDUSP_n_ppERK_n + 220.0*p52*pDUSP_n*ppERK_n;
    w[47] = 220.0*p53*pDUSP_n_ppERK_n;
    w[48] = -220.0*m54*pDUSP_n_pERK_n + 220.0*p54*pDUSP_n*pERK_n;
    w[49] = 220.0*p55*pDUSP_n_pERK_n;
    w[50] = 220.0*ERK_n*p56*pDUSP_n - 220.0*m56*pDUSP_n_ERK_n;
    w[51] = 220.0*DUSP_n*p47*ppERK_n - 220.0*DUSP_n_ppERK_n*m47;
    w[52] = 220.0*DUSP_n_ppERK_n*p48;
    w[53] = 220.0*DUSP_n*p49*pERK_n - 220.0*DUSP_n_pERK_n*m49;
    w[54] = 220.0*DUSP_n_pERK_n*p50;
    w[55] = 220.0*DUSP_n*ERK_n*p51 - 220.0*DUSP_n_ERK_n*m51;
    w[56] = 220.0*pow(FOSn_2, n57)*V57/(pow(FOSn_2, n57) + pow(K57, n57));
    w[57] = PreFmRNA*Vn*p58;
    w[58] = 940.0*FmRNA*p59;
    w[59] = 940.0*FmRNA*p60;
    w[60] = 940.0*F*p61;
    w[61] = F*KimF*Vc - Fn*KexF*Vn;
    w[62] = 940.0*Fn*p63;
    w[63] = 940.0*A1*EGF*V101/(A1 + K101);
    w[64] = 940.0*A1_2*V102/(A1_2 + K102);
    w[65] = 940.0*A2*HRG*V103/(A2 + K103);
    w[66] = 940.0*A2_2*V104/(A2_2 + K104);
    w[67] = 940.0*EGF*RsD*V105/(K105 + RsD);
    w[68] = 940.0*HRG*RsD*V106/(K106 + RsD);
    w[69] = 940.0*A1_2*RsT*V107/(K107 + RsT);
    w[70] = 940.0*A2_2*RsT*V108/(K108 + RsT);
    w[71] = 940.0*A3*HRG*V109/(A3 + K109);
    w[72] = 940.0*A3_2*V110/(A3_2 + K110);
    w[73] = 940.0*HRG*Kin*V111/(K111 + Kin);
    w[74] = 940.0*A3_2*Kin_2*V112/(K112 + Kin_2);
    w[75] = 940.0*MEK*RsT*V113/(K113 + MEK);
    w[76] = 940.0*Kin_2*MEK*V114/(K114 + MEK);
    w[77] = 940.0*V115*pMEK/(K115 + pMEK);
}