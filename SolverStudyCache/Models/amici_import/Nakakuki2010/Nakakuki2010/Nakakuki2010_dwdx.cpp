#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Nakakuki2010(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 940.0*A1*V101/(A1 + K101);
    dwdx[1] = 940.0*RsD*V105/(K105 + RsD);
    dwdx[2] = 940.0*A2*V103/(A2 + K103);
    dwdx[3] = 940.0*RsD*V106/(K106 + RsD);
    dwdx[4] = 940.0*A3*V109/(A3 + K109);
    dwdx[5] = 940.0*Kin*V111/(K111 + Kin);
    dwdx[6] = -940.0*A1*EGF*V101/pow(A1 + K101, 2) + 940.0*EGF*V101/(A1 + K101);
    dwdx[7] = -940.0*A1_2*V102/pow(A1_2 + K102, 2) + 940.0*V102/(A1_2 + K102);
    dwdx[8] = 940.0*RsT*V107/(K107 + RsT);
    dwdx[9] = -940.0*A2*HRG*V103/pow(A2 + K103, 2) + 940.0*HRG*V103/(A2 + K103);
    dwdx[10] = -940.0*A2_2*V104/pow(A2_2 + K104, 2) + 940.0*V104/(A2_2 + K104);
    dwdx[11] = 940.0*RsT*V108/(K108 + RsT);
    dwdx[12] = -940.0*A3*HRG*V109/pow(A3 + K109, 2) + 940.0*HRG*V109/(A3 + K109);
    dwdx[13] = -940.0*A3_2*V110/pow(A3_2 + K110, 2) + 940.0*V110/(A3_2 + K110);
    dwdx[14] = 940.0*Kin_2*V112/(K112 + Kin_2);
    dwdx[15] = 940.0*p12;
    dwdx[16] = 940.0*p13;
    dwdx[17] = -940.0*ERK_c*Fct*V1*pMEK/pow(ERK_c + K1*(1 + pERK_c/K2), 2) + 940.0*Fct*V1*pMEK/(ERK_c + K1*(1 + pERK_c/K2));
    dwdx[18] = -940.0*Fct*K2*V2*pERK_c*pMEK/(K1*pow(K2*(ERK_c/K1 + 1) + pERK_c, 2));
    dwdx[19] = KimERK*Vc;
    dwdx[20] = -940.0*ERK_c*Fct*K1*V1*pMEK/(K2*pow(ERK_c + K1*(1 + pERK_c/K2), 2));
    dwdx[21] = -940.0*Fct*V2*pERK_c*pMEK/pow(K2*(ERK_c/K1 + 1) + pERK_c, 2) + 940.0*Fct*V2*pMEK/(K2*(ERK_c/K1 + 1) + pERK_c);
    dwdx[22] = -940.0*V3*pERK_c/pow(K3*(1 + ppERK_c/K4) + pERK_c, 2) + 940.0*V3/(K3*(1 + ppERK_c/K4) + pERK_c);
    dwdx[23] = -940.0*K4*V4*ppERK_c/(K3*pow(K4*(1 + pERK_c/K3) + ppERK_c, 2));
    dwdx[24] = KimERKP*Vc;
    dwdx[25] = -940.0*K3*V3*pERK_c/(K4*pow(K3*(1 + ppERK_c/K4) + pERK_c, 2));
    dwdx[26] = -940.0*V4*ppERK_c/pow(K4*(1 + pERK_c/K3) + ppERK_c, 2) + 940.0*V4/(K4*(1 + pERK_c/K3) + ppERK_c);
    dwdx[27] = KimERKPP*Vc;
    dwdx[28] = 940.0*DUSP_c*V14/(DUSP_c + K14);
    dwdx[29] = 940.0*RSK_c*V24/(K24 + RSK_c);
    dwdx[30] = 940.0*V35*c_FOS_c/(K35 + c_FOS_c);
    dwdx[31] = 940.0*p61;
    dwdx[32] = KimF*Vc;
    dwdx[33] = -940.0*V35*c_FOS_c*ppERK_c/pow(K35 + c_FOS_c, 2) + 940.0*V35*ppERK_c/(K35 + c_FOS_c);
    dwdx[34] = -940.0*V36*c_FOS_c*pRSK_c/pow(K36 + c_FOS_c, 2) + 940.0*V36*pRSK_c/(K36 + c_FOS_c);
    dwdx[35] = 940.0*p38;
    dwdx[36] = KimFOS*Vc;
    dwdx[37] = -940.0*V37*pc_FOS_c/pow(K37 + pc_FOS_c, 2) + 940.0*V37/(K37 + pc_FOS_c);
    dwdx[38] = 940.0*p39;
    dwdx[39] = KimFOSP*Vc;
    dwdx[40] = 940.0*p33;
    dwdx[41] = 940.0*p34;
    dwdx[42] = 940.0*p59;
    dwdx[43] = 940.0*p60;
    dwdx[44] = -940.0*HRG*Kin*V111/pow(K111 + Kin, 2) + 940.0*HRG*V111/(K111 + Kin);
    dwdx[45] = -940.0*A3_2*Kin_2*V112/pow(K112 + Kin_2, 2) + 940.0*A3_2*V112/(K112 + Kin_2);
    dwdx[46] = 940.0*MEK*V114/(K114 + MEK);
    dwdx[47] = 940.0*ERK_c*Fct*V1/(ERK_c + K1*(1 + pERK_c/K2));
    dwdx[48] = 940.0*Fct*V2*pERK_c/(K2*(ERK_c/K1 + 1) + pERK_c);
    dwdx[49] = -940.0*V115*pMEK/pow(K115 + pMEK, 2) + 940.0*V115/(K115 + pMEK);
    dwdx[50] = -940.0*MEK*RsT*V113/pow(K113 + MEK, 2) + 940.0*RsT*V113/(K113 + MEK);
    dwdx[51] = -940.0*Kin_2*MEK*V114/pow(K114 + MEK, 2) + 940.0*Kin_2*V114/(K114 + MEK);
    dwdx[52] = -940.0*DUSP_c*V14*ppERK_c/pow(DUSP_c + K14, 2) + 940.0*V14*ppERK_c/(DUSP_c + K14);
    dwdx[53] = 940.0*p16;
    dwdx[54] = KimDUSP*Vc;
    dwdx[55] = -940.0*V15*pDUSP_c/pow(K15 + pDUSP_c, 2) + 940.0*V15/(K15 + pDUSP_c);
    dwdx[56] = 940.0*p17;
    dwdx[57] = KimDUSPP*Vc;
    dwdx[58] = -940.0*RSK_c*V24*ppERK_c/pow(K24 + RSK_c, 2) + 940.0*V24*ppERK_c/(K24 + RSK_c);
    dwdx[59] = -940.0*V25*pRSK_c/pow(K25 + pRSK_c, 2) + 940.0*V25/(K25 + pRSK_c);
    dwdx[60] = KimRSKP*Vc;
    dwdx[61] = 940.0*V36*c_FOS_c/(K36 + c_FOS_c);
    dwdx[62] = -940.0*EGF*RsD*V105/pow(K105 + RsD, 2) + 940.0*EGF*V105/(K105 + RsD);
    dwdx[63] = -940.0*HRG*RsD*V106/pow(K106 + RsD, 2) + 940.0*HRG*V106/(K106 + RsD);
    dwdx[64] = -940.0*A1_2*RsT*V107/pow(K107 + RsT, 2) + 940.0*A1_2*V107/(K107 + RsT);
    dwdx[65] = -940.0*A2_2*RsT*V108/pow(K108 + RsT, 2) + 940.0*A2_2*V108/(K108 + RsT);
    dwdx[66] = 940.0*MEK*V113/(K113 + MEK);
    dwdx[67] = -220.0*CREB_n*V27*pRSK_n/pow(CREB_n + K27, 2) + 220.0*V27*pRSK_n/(CREB_n + K27);
    dwdx[68] = -220.0*V28*pCREB_n/pow(K28 + pCREB_n, 2) + 220.0*V28/(K28 + pCREB_n);
    dwdx[69] = -220.0*V31*n31*pow(pCREB_n*pElk1_n, 2*n31)/(pCREB_n*pow(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31), 2)) + 220.0*V31*n31*pow(pCREB_n*pElk1_n, n31)/(pCREB_n*(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31)));
    dwdx[70] = -KexERK*Vn;
    dwdx[71] = 220.0*p56*pDUSP_n;
    dwdx[72] = 220.0*DUSP_n*p51;
    dwdx[73] = -220.0*V5*pERK_n/pow(K5*(1 + ppERK_n/K6) + pERK_n, 2) + 220.0*V5/(K5*(1 + ppERK_n/K6) + pERK_n);
    dwdx[74] = -220.0*K6*V6*ppERK_n/(K5*pow(K6*(1 + pERK_n/K5) + ppERK_n, 2));
    dwdx[75] = -KexERKP*Vn;
    dwdx[76] = 220.0*p54*pDUSP_n;
    dwdx[77] = 220.0*DUSP_n*p49;
    dwdx[78] = -220.0*K5*V5*pERK_n/(K6*pow(K5*(1 + ppERK_n/K6) + pERK_n, 2));
    dwdx[79] = -220.0*V6*ppERK_n/pow(K6*(1 + pERK_n/K5) + ppERK_n, 2) + 220.0*V6/(K6*(1 + pERK_n/K5) + ppERK_n);
    dwdx[80] = -KexERKPP*Vn;
    dwdx[81] = -220.0*V10*n10*pow(ppERK_n, 2*n10)/(ppERK_n*pow(pow(K10, n10) + pow(ppERK_n, n10), 2)) + 220.0*V10*n10*pow(ppERK_n, n10)/(ppERK_n*(pow(K10, n10) + pow(ppERK_n, n10)));
    dwdx[82] = 220.0*DUSP_n*V20/(DUSP_n + K20);
    dwdx[83] = 220.0*Elk1_n*V29/(Elk1_n + K29);
    dwdx[84] = 220.0*FOSn*V42/(FOSn + K42);
    dwdx[85] = 220.0*p52*pDUSP_n;
    dwdx[86] = 220.0*DUSP_n*p47;
    dwdx[87] = -220.0*Elk1_n*V29*ppERK_n/pow(Elk1_n + K29, 2) + 220.0*V29*ppERK_n/(Elk1_n + K29);
    dwdx[88] = -220.0*V30*pElk1_n/pow(K30 + pElk1_n, 2) + 220.0*V30/(K30 + pElk1_n);
    dwdx[89] = -220.0*V31*n31*pow(pCREB_n*pElk1_n, 2*n31)/(pElk1_n*pow(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31), 2)) + 220.0*V31*n31*pow(pCREB_n*pElk1_n, n31)/(pElk1_n*(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31)));
    dwdx[90] = -KexFOS*Vn;
    dwdx[91] = -220.0*FOSn*V42*ppERK_n/pow(FOSn + K42, 2) + 220.0*V42*ppERK_n/(FOSn + K42);
    dwdx[92] = -220.0*FOSn*V43*pRSK_n/pow(FOSn + K43, 2) + 220.0*V43*pRSK_n/(FOSn + K43);
    dwdx[93] = 220.0*p45;
    dwdx[94] = -KexFOSP*Vn;
    dwdx[95] = -220.0*FOSn_2*V44/pow(FOSn_2 + K44, 2) + 220.0*V44/(FOSn_2 + K44);
    dwdx[96] = 220.0*p46;
    dwdx[97] = -220.0*pow(FOSn_2, 2*n57)*V57*n57/(FOSn_2*pow(pow(FOSn_2, n57) + pow(K57, n57), 2)) + 220.0*pow(FOSn_2, n57)*V57*n57/(FOSn_2*(pow(FOSn_2, n57) + pow(K57, n57)));
    dwdx[98] = -220.0*V31*nF31*pow(Fn/KF31, nF31)*pow(pCREB_n*pElk1_n, n31)/(Fn*pow(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31), 2));
    dwdx[99] = -KexF*Vn;
    dwdx[100] = 940.0*p63;
    dwdx[101] = -KexDUSP*Vn;
    dwdx[102] = -220.0*DUSP_n*V20*ppERK_n/pow(DUSP_n + K20, 2) + 220.0*V20*ppERK_n/(DUSP_n + K20);
    dwdx[103] = 220.0*p22;
    dwdx[104] = 220.0*p47*ppERK_n;
    dwdx[105] = 220.0*p49*pERK_n;
    dwdx[106] = 220.0*ERK_n*p51;
    dwdx[107] = -KexDUSPP*Vn;
    dwdx[108] = -220.0*V21*pDUSP_n/pow(K21 + pDUSP_n, 2) + 220.0*V21/(K21 + pDUSP_n);
    dwdx[109] = 220.0*p23;
    dwdx[110] = 220.0*p52*ppERK_n;
    dwdx[111] = 220.0*p54*pERK_n;
    dwdx[112] = 220.0*ERK_n*p56;
    dwdx[113] = -220.0*m56;
    dwdx[114] = -220.0*m54;
    dwdx[115] = 220.0*p55;
    dwdx[116] = -220.0*m52;
    dwdx[117] = 220.0*p53;
    dwdx[118] = -220.0*m51;
    dwdx[119] = -220.0*m49;
    dwdx[120] = 220.0*p50;
    dwdx[121] = -220.0*m47;
    dwdx[122] = 220.0*p48;
    dwdx[123] = Vn*p11;
    dwdx[124] = Vn*p32;
    dwdx[125] = Vn*p58;
    dwdx[126] = -KexRSKP*Vn;
    dwdx[127] = 220.0*CREB_n*V27/(CREB_n + K27);
    dwdx[128] = 220.0*FOSn*V43/(FOSn + K43);
}