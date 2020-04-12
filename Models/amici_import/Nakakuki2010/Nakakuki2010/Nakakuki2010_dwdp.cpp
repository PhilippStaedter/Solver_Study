#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Nakakuki2010(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = 940.0*Fct*pERK_c*pMEK/(K2*(ERK_c/K1 + 1) + pERK_c);
            break;
        case 1:
            dwdp[2] = 940.0*pERK_c/(K3*(1 + ppERK_c/K4) + pERK_c);
            break;
        case 2:
            dwdp[2] = 940.0*V3*pERK_c*(-1 - ppERK_c/K4)/pow(K3*(1 + ppERK_c/K4) + pERK_c, 2);
            dwdp[3] = 940.0*K4*V4*pERK_c*ppERK_c/(pow(K3, 2)*pow(K4*(1 + pERK_c/K3) + ppERK_c, 2));
            break;
        case 3:
            dwdp[3] = 940.0*ppERK_c/(K4*(1 + pERK_c/K3) + ppERK_c);
            break;
        case 4:
            dwdp[2] = 940.0*K3*V3*pERK_c*ppERK_c/(pow(K4, 2)*pow(K3*(1 + ppERK_c/K4) + pERK_c, 2));
            dwdp[3] = 940.0*V4*ppERK_c*(-1 - pERK_c/K3)/pow(K4*(1 + pERK_c/K3) + ppERK_c, 2);
            break;
        case 5:
            dwdp[4] = 220.0*pERK_n/(K5*(1 + ppERK_n/K6) + pERK_n);
            break;
        case 6:
            dwdp[4] = 220.0*V5*pERK_n*(-1 - ppERK_n/K6)/pow(K5*(1 + ppERK_n/K6) + pERK_n, 2);
            dwdp[5] = 220.0*K6*V6*pERK_n*ppERK_n/(pow(K5, 2)*pow(K6*(1 + pERK_n/K5) + ppERK_n, 2));
            break;
        case 7:
            dwdp[5] = 220.0*ppERK_n/(K6*(1 + pERK_n/K5) + ppERK_n);
            break;
        case 8:
            dwdp[4] = 220.0*K5*V5*pERK_n*ppERK_n/(pow(K6, 2)*pow(K5*(1 + ppERK_n/K6) + pERK_n, 2));
            dwdp[5] = 220.0*V6*ppERK_n*(-1 - pERK_n/K5)/pow(K6*(1 + pERK_n/K5) + ppERK_n, 2);
            break;
        case 9:
            dwdp[6] = ERK_c*Vc;
            break;
        case 10:
            dwdp[6] = -ERK_n*Vn;
            break;
        case 11:
            dwdp[7] = Vc*pERK_c;
            break;
        case 12:
            dwdp[7] = -Vn*pERK_n;
            break;
        case 13:
            dwdp[8] = Vc*ppERK_c;
            break;
        case 14:
            dwdp[8] = -Vn*ppERK_n;
            break;
        case 15:
            dwdp[9] = 220.0*pow(ppERK_n, n10)/(pow(K10, n10) + pow(ppERK_n, n10));
            break;
        case 16:
            dwdp[9] = -220.0*pow(K10, n10)*V10*n10*pow(ppERK_n, n10)/(K10*pow(pow(K10, n10) + pow(ppERK_n, n10), 2));
            break;
        case 17:
            dwdp[9] = 220.0*V10*pow(ppERK_n, n10)*log(ppERK_n)/(pow(K10, n10) + pow(ppERK_n, n10)) + 220.0*V10*pow(ppERK_n, n10)*(-pow(K10, n10)*log(K10) - pow(ppERK_n, n10)*log(ppERK_n))/pow(pow(K10, n10) + pow(ppERK_n, n10), 2);
            break;
        case 18:
            dwdp[10] = PreDUSPmRNA*Vn;
            break;
        case 19:
            dwdp[11] = 940.0*DUSPmRNA;
            break;
        case 20:
            dwdp[12] = 940.0*DUSPmRNA;
            break;
        case 21:
            dwdp[13] = 940.0*DUSP_c*ppERK_c/(DUSP_c + K14);
            break;
        case 22:
            dwdp[13] = -940.0*DUSP_c*V14*ppERK_c/pow(DUSP_c + K14, 2);
            break;
        case 23:
            dwdp[14] = 940.0*pDUSP_c/(K15 + pDUSP_c);
            break;
        case 24:
            dwdp[14] = -940.0*V15*pDUSP_c/pow(K15 + pDUSP_c, 2);
            break;
        case 25:
            dwdp[15] = 940.0*DUSP_c;
            break;
        case 26:
            dwdp[16] = 940.0*pDUSP_c;
            break;
        case 27:
            dwdp[17] = DUSP_c*Vc;
            break;
        case 28:
            dwdp[17] = -DUSP_n*Vn;
            break;
        case 29:
            dwdp[18] = Vc*pDUSP_c;
            break;
        case 30:
            dwdp[18] = -Vn*pDUSP_n;
            break;
        case 31:
            dwdp[19] = 220.0*DUSP_n*ppERK_n/(DUSP_n + K20);
            break;
        case 32:
            dwdp[19] = -220.0*DUSP_n*V20*ppERK_n/pow(DUSP_n + K20, 2);
            break;
        case 33:
            dwdp[20] = 220.0*pDUSP_n/(K21 + pDUSP_n);
            break;
        case 34:
            dwdp[20] = -220.0*V21*pDUSP_n/pow(K21 + pDUSP_n, 2);
            break;
        case 35:
            dwdp[21] = 220.0*DUSP_n;
            break;
        case 36:
            dwdp[22] = 220.0*pDUSP_n;
            break;
        case 37:
            dwdp[23] = 940.0*RSK_c*ppERK_c/(K24 + RSK_c);
            break;
        case 38:
            dwdp[23] = -940.0*RSK_c*V24*ppERK_c/pow(K24 + RSK_c, 2);
            break;
        case 39:
            dwdp[0] = 940.0*ERK_c*Fct*pMEK/(ERK_c + K1*(1 + pERK_c/K2));
            break;
        case 40:
            dwdp[0] = 940.0*ERK_c*Fct*V1*pMEK*(-1 - pERK_c/K2)/pow(ERK_c + K1*(1 + pERK_c/K2), 2);
            dwdp[1] = 940.0*ERK_c*Fct*K2*V2*pERK_c*pMEK/(pow(K1, 2)*pow(K2*(ERK_c/K1 + 1) + pERK_c, 2));
            break;
        case 41:
            dwdp[24] = 940.0*pRSK_c/(K25 + pRSK_c);
            break;
        case 42:
            dwdp[24] = -940.0*V25*pRSK_c/pow(K25 + pRSK_c, 2);
            break;
        case 43:
            dwdp[25] = Vc*pRSK_c;
            break;
        case 44:
            dwdp[25] = -Vn*pRSK_n;
            break;
        case 45:
            dwdp[26] = 220.0*CREB_n*pRSK_n/(CREB_n + K27);
            break;
        case 46:
            dwdp[26] = -220.0*CREB_n*V27*pRSK_n/pow(CREB_n + K27, 2);
            break;
        case 47:
            dwdp[27] = 220.0*pCREB_n/(K28 + pCREB_n);
            break;
        case 48:
            dwdp[27] = -220.0*V28*pCREB_n/pow(K28 + pCREB_n, 2);
            break;
        case 49:
            dwdp[28] = 220.0*Elk1_n*ppERK_n/(Elk1_n + K29);
            break;
        case 50:
            dwdp[28] = -220.0*Elk1_n*V29*ppERK_n/pow(Elk1_n + K29, 2);
            break;
        case 51:
            dwdp[29] = 220.0*pElk1_n/(K30 + pElk1_n);
            break;
        case 52:
            dwdp[29] = -220.0*V30*pElk1_n/pow(K30 + pElk1_n, 2);
            break;
        case 53:
            dwdp[30] = 220.0*pow(pCREB_n*pElk1_n, n31)/(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31));
            break;
        case 54:
            dwdp[30] = -220.0*pow(K31, n31)*V31*n31*pow(pCREB_n*pElk1_n, n31)/(K31*pow(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31), 2));
            break;
        case 55:
            dwdp[30] = 220.0*V31*pow(pCREB_n*pElk1_n, n31)*(-pow(K31, n31)*log(K31) - pow(pCREB_n*pElk1_n, n31)*log(pCREB_n*pElk1_n))/pow(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31), 2) + 220.0*V31*pow(pCREB_n*pElk1_n, n31)*log(pCREB_n*pElk1_n)/(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31));
            break;
        case 56:
            dwdp[31] = PreFOSmRNA*Vn;
            break;
        case 57:
            dwdp[32] = 940.0*c_FOSmRNA;
            break;
        case 58:
            dwdp[33] = 940.0*c_FOSmRNA;
            break;
        case 59:
            dwdp[34] = 940.0*c_FOS_c*ppERK_c/(K35 + c_FOS_c);
            break;
        case 60:
            dwdp[34] = -940.0*V35*c_FOS_c*ppERK_c/pow(K35 + c_FOS_c, 2);
            break;
        case 61:
            dwdp[35] = 940.0*c_FOS_c*pRSK_c/(K36 + c_FOS_c);
            break;
        case 62:
            dwdp[35] = -940.0*V36*c_FOS_c*pRSK_c/pow(K36 + c_FOS_c, 2);
            break;
        case 63:
            dwdp[36] = 940.0*pc_FOS_c/(K37 + pc_FOS_c);
            break;
        case 64:
            dwdp[36] = -940.0*V37*pc_FOS_c/pow(K37 + pc_FOS_c, 2);
            break;
        case 65:
            dwdp[37] = 940.0*c_FOS_c;
            break;
        case 66:
            dwdp[38] = 940.0*pc_FOS_c;
            break;
        case 67:
            dwdp[39] = Vc*c_FOS_c;
            break;
        case 68:
            dwdp[39] = -FOSn*Vn;
            break;
        case 69:
            dwdp[40] = Vc*pc_FOS_c;
            break;
        case 70:
            dwdp[40] = -FOSn_2*Vn;
            break;
        case 71:
            dwdp[41] = 220.0*FOSn*ppERK_n/(FOSn + K42);
            break;
        case 72:
            dwdp[41] = -220.0*FOSn*V42*ppERK_n/pow(FOSn + K42, 2);
            break;
        case 73:
            dwdp[42] = 220.0*FOSn*pRSK_n/(FOSn + K43);
            break;
        case 74:
            dwdp[42] = -220.0*FOSn*V43*pRSK_n/pow(FOSn + K43, 2);
            break;
        case 75:
            dwdp[43] = 220.0*FOSn_2/(FOSn_2 + K44);
            break;
        case 76:
            dwdp[43] = -220.0*FOSn_2*V44/pow(FOSn_2 + K44, 2);
            break;
        case 77:
            dwdp[44] = 220.0*FOSn;
            break;
        case 78:
            dwdp[45] = 220.0*FOSn_2;
            break;
        case 79:
            dwdp[51] = 220.0*DUSP_n*ppERK_n;
            break;
        case 80:
            dwdp[51] = -220.0*DUSP_n_ppERK_n;
            break;
        case 81:
            dwdp[52] = 220.0*DUSP_n_ppERK_n;
            break;
        case 82:
            dwdp[53] = 220.0*DUSP_n*pERK_n;
            break;
        case 83:
            dwdp[53] = -220.0*DUSP_n_pERK_n;
            break;
        case 84:
            dwdp[54] = 220.0*DUSP_n_pERK_n;
            break;
        case 85:
            dwdp[55] = 220.0*DUSP_n*ERK_n;
            break;
        case 86:
            dwdp[55] = -220.0*DUSP_n_ERK_n;
            break;
        case 87:
            dwdp[0] = 940.0*ERK_c*V1*pMEK/(ERK_c + K1*(1 + pERK_c/K2));
            dwdp[1] = 940.0*V2*pERK_c*pMEK/(K2*(ERK_c/K1 + 1) + pERK_c);
            break;
        case 88:
            dwdp[46] = 220.0*pDUSP_n*ppERK_n;
            break;
        case 89:
            dwdp[46] = -220.0*pDUSP_n_ppERK_n;
            break;
        case 90:
            dwdp[47] = 220.0*pDUSP_n_ppERK_n;
            break;
        case 91:
            dwdp[48] = 220.0*pDUSP_n*pERK_n;
            break;
        case 92:
            dwdp[48] = -220.0*pDUSP_n_pERK_n;
            break;
        case 93:
            dwdp[49] = 220.0*pDUSP_n_pERK_n;
            break;
        case 94:
            dwdp[50] = 220.0*ERK_n*pDUSP_n;
            break;
        case 95:
            dwdp[50] = -220.0*pDUSP_n_ERK_n;
            break;
        case 96:
            dwdp[56] = 220.0*pow(FOSn_2, n57)/(pow(FOSn_2, n57) + pow(K57, n57));
            break;
        case 97:
            dwdp[56] = -220.0*pow(FOSn_2, n57)*pow(K57, n57)*V57*n57/(K57*pow(pow(FOSn_2, n57) + pow(K57, n57), 2));
            break;
        case 98:
            dwdp[56] = 220.0*pow(FOSn_2, n57)*V57*log(FOSn_2)/(pow(FOSn_2, n57) + pow(K57, n57)) + 220.0*pow(FOSn_2, n57)*V57*(-pow(FOSn_2, n57)*log(FOSn_2) - pow(K57, n57)*log(K57))/pow(pow(FOSn_2, n57) + pow(K57, n57), 2);
            break;
        case 99:
            dwdp[57] = PreFmRNA*Vn;
            break;
        case 100:
            dwdp[58] = 940.0*FmRNA;
            break;
        case 101:
            dwdp[59] = 940.0*FmRNA;
            break;
        case 102:
            dwdp[60] = 940.0*F;
            break;
        case 103:
            dwdp[61] = F*Vc;
            break;
        case 104:
            dwdp[61] = -Fn*Vn;
            break;
        case 105:
            dwdp[62] = 940.0*Fn;
            break;
        case 106:
            dwdp[30] = 220.0*V31*nF31*pow(Fn/KF31, nF31)*pow(pCREB_n*pElk1_n, n31)/(KF31*pow(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31), 2));
            break;
        case 107:
            dwdp[30] = -220.0*V31*pow(Fn/KF31, nF31)*pow(pCREB_n*pElk1_n, n31)*log(Fn/KF31)/pow(pow(K31, n31) + pow(Fn/KF31, nF31) + pow(pCREB_n*pElk1_n, n31), 2);
            break;
        case 108:
            dwdp[0] = 940.0*ERK_c*Fct*K1*V1*pERK_c*pMEK/(pow(K2, 2)*pow(ERK_c + K1*(1 + pERK_c/K2), 2));
            dwdp[1] = 940.0*Fct*V2*pERK_c*pMEK*(-ERK_c/K1 - 1)/pow(K2*(ERK_c/K1 + 1) + pERK_c, 2);
            break;
        case 109:
            dwdp[6] = -ERK_n*KexERK;
            dwdp[7] = -KexERKP*pERK_n;
            dwdp[8] = -KexERKPP*ppERK_n;
            dwdp[10] = PreDUSPmRNA*p11;
            dwdp[17] = -DUSP_n*KexDUSP;
            dwdp[18] = -KexDUSPP*pDUSP_n;
            dwdp[25] = -KexRSKP*pRSK_n;
            dwdp[31] = PreFOSmRNA*p32;
            dwdp[39] = -FOSn*KexFOS;
            dwdp[40] = -FOSn_2*KexFOSP;
            dwdp[57] = PreFmRNA*p58;
            dwdp[61] = -Fn*KexF;
            break;
        case 110:
            dwdp[6] = ERK_c*KimERK;
            dwdp[7] = KimERKP*pERK_c;
            dwdp[8] = KimERKPP*ppERK_c;
            dwdp[17] = DUSP_c*KimDUSP;
            dwdp[18] = KimDUSPP*pDUSP_c;
            dwdp[25] = KimRSKP*pRSK_c;
            dwdp[39] = KimFOS*c_FOS_c;
            dwdp[40] = KimFOSP*pc_FOS_c;
            dwdp[61] = F*KimF;
            break;
        case 111:
            dwdp[63] = 940.0*A1*EGF/(A1 + K101);
            break;
        case 112:
            dwdp[63] = -940.0*A1*EGF*V101/pow(A1 + K101, 2);
            break;
        case 113:
            dwdp[64] = 940.0*A1_2/(A1_2 + K102);
            break;
        case 114:
            dwdp[64] = -940.0*A1_2*V102/pow(A1_2 + K102, 2);
            break;
        case 115:
            dwdp[65] = 940.0*A2*HRG/(A2 + K103);
            break;
        case 116:
            dwdp[65] = -940.0*A2*HRG*V103/pow(A2 + K103, 2);
            break;
        case 117:
            dwdp[66] = 940.0*A2_2/(A2_2 + K104);
            break;
        case 118:
            dwdp[66] = -940.0*A2_2*V104/pow(A2_2 + K104, 2);
            break;
        case 119:
            dwdp[67] = 940.0*EGF*RsD/(K105 + RsD);
            break;
        case 120:
            dwdp[67] = -940.0*EGF*RsD*V105/pow(K105 + RsD, 2);
            break;
        case 121:
            dwdp[68] = 940.0*HRG*RsD/(K106 + RsD);
            break;
        case 122:
            dwdp[68] = -940.0*HRG*RsD*V106/pow(K106 + RsD, 2);
            break;
        case 123:
            dwdp[69] = 940.0*A1_2*RsT/(K107 + RsT);
            break;
        case 124:
            dwdp[69] = -940.0*A1_2*RsT*V107/pow(K107 + RsT, 2);
            break;
        case 125:
            dwdp[70] = 940.0*A2_2*RsT/(K108 + RsT);
            break;
        case 126:
            dwdp[70] = -940.0*A2_2*RsT*V108/pow(K108 + RsT, 2);
            break;
        case 127:
            dwdp[71] = 940.0*A3*HRG/(A3 + K109);
            break;
        case 128:
            dwdp[71] = -940.0*A3*HRG*V109/pow(A3 + K109, 2);
            break;
        case 129:
            dwdp[72] = 940.0*A3_2/(A3_2 + K110);
            break;
        case 130:
            dwdp[72] = -940.0*A3_2*V110/pow(A3_2 + K110, 2);
            break;
        case 131:
            dwdp[73] = 940.0*HRG*Kin/(K111 + Kin);
            break;
        case 132:
            dwdp[73] = -940.0*HRG*Kin*V111/pow(K111 + Kin, 2);
            break;
        case 133:
            dwdp[74] = 940.0*A3_2*Kin_2/(K112 + Kin_2);
            break;
        case 134:
            dwdp[74] = -940.0*A3_2*Kin_2*V112/pow(K112 + Kin_2, 2);
            break;
        case 135:
            dwdp[75] = 940.0*MEK*RsT/(K113 + MEK);
            break;
        case 136:
            dwdp[75] = -940.0*MEK*RsT*V113/pow(K113 + MEK, 2);
            break;
        case 137:
            dwdp[76] = 940.0*Kin_2*MEK/(K114 + MEK);
            break;
        case 138:
            dwdp[76] = -940.0*Kin_2*MEK*V114/pow(K114 + MEK, 2);
            break;
        case 139:
            dwdp[77] = 940.0*pMEK/(K115 + pMEK);
            break;
        case 140:
            dwdp[77] = -940.0*V115*pMEK/pow(K115 + pMEK, 2);
            break;
    }
}