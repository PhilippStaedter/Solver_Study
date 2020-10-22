#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Proctor2013a(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[42] = 1.0*MMP3*proMMP13;
            break;
        case 1:
            dwdp[29] = 1.0*Matriptase*proMMP1;
            break;
        case 2:
            dwdp[30] = 1.0*MMP3*proMMP1;
            break;
        case 3:
            dwdp[36] = 1.0*Matriptase*proMMP3;
            break;
        case 4:
            dwdp[16] = 1.0*cFos_cJun*ksyncJunmRNA;
            dwdp[25] = 1.0*cFos_cJun*ksynMMP1mRNA;
            dwdp[32] = 1.0*cFos_cJun*ksynMMP3mRNA;
            dwdp[38] = 1.0*cFos_cJun*ksynMMP13mRNA;
            dwdp[47] = 1.0*cFos_cJun*ksynADAMTS4mRNA;
            dwdp[52] = 1.0*cFos_cJun*ksynPP4;
            dwdp[54] = 1.0*cFos_cJun*ksynDUSP16;
            dwdp[62] = 1.0*cFos_cJun*ksyncFosmRNA;
            dwdp[68] = 1.0*cFos_cJun*ksynMKP1;
            dwdp[109] = 1.0*cFos_cJun*ksynMatriptase;
            dwdp[111] = 1.0*cFos_cJun*ksynSP1;
            dwdp[117] = 1.0*TIMP1_DNA*cFos_cJun*ksynTIMP1mRNA;
            dwdp[119] = 1.0*cFos_cJun*ksynTIMP3mRNA;
            dwdp[120] = 1.0*STAT3_P_nuc*ksynTIMP3mRNAStat3;
            break;
        case 5:
            dwdp[72] = 1.0*cFos_P*cJun_P;
            break;
        case 6:
            dwdp[0] = 1.0*IL1*IL1R;
            break;
        case 7:
            dwdp[2] = 1.0*IL1*IL1Ra;
            break;
        case 8:
            dwdp[5] = 1.0*IL1_IL1R*IRAK2;
            break;
        case 9:
            dwdp[74] = 1.0*OSM*OSMR;
            break;
        case 10:
            dwdp[76] = 1.0*OSM*OSMRa;
            break;
        case 11:
            dwdp[95] = 1.0*OSMR*SOCS3;
            break;
        case 12:
            dwdp[113] = 1.0*SP1*TIMP1_DNA;
            break;
        case 13:
            dwdp[7] = 1.0*IL1_IL1R_IRAK2*TRAF6;
            break;
        case 14:
            dwdp[86] = 1.0*STAT3_P_cyt;
            break;
        case 15:
            dwdp[15] = 1.0*cJun_dimer;
            break;
        case 16:
            dwdp[51] = 1.0*ADAMTS4;
            break;
        case 17:
            dwdp[50] = 1.0*ADAMTS4_mRNA;
            break;
        case 18:
            dwdp[108] = 1.0*ADAMTS4*Aggrecan_Collagen2;
            break;
        case 19:
            dwdp[65] = 1.0*cFos;
            break;
        case 20:
            dwdp[63] = 1.0*cFos_mRNA;
            break;
        case 21:
            dwdp[21] = 1.0*cJun;
            break;
        case 22:
            dwdp[19] = 1.0*cJun_mRNA;
            break;
        case 23:
            dwdp[106] = 1.0*Collagen2*MMP1;
            break;
        case 24:
            dwdp[107] = 1.0*Collagen2*MMP13;
            break;
        case 25:
            dwdp[57] = 1.0*DUSP16;
            break;
        case 26:
            dwdp[4] = 1.0*IL1;
            break;
        case 27:
            dwdp[110] = 1.0*Matriptase;
            break;
        case 28:
            dwdp[70] = 1.0*MKP1;
            break;
        case 29:
            dwdp[31] = 1.0*MMP1;
            break;
        case 30:
            dwdp[43] = 1.0*MMP13;
            break;
        case 31:
            dwdp[41] = 1.0*MMP13_mRNA;
            break;
        case 32:
            dwdp[28] = 1.0*MMP1_mRNA;
            break;
        case 33:
            dwdp[37] = 1.0*MMP3;
            break;
        case 34:
            dwdp[35] = 1.0*MMP3_mRNA;
            break;
        case 35:
            dwdp[97] = 1.0*OSM;
            break;
        case 36:
            dwdp[56] = 1.0*PP4;
            break;
        case 37:
            dwdp[90] = 1.0*PTPRT;
            break;
        case 38:
            dwdp[94] = 1.0*SOCS3;
            break;
        case 39:
            dwdp[92] = 1.0*SOCS3_mRNA;
            break;
        case 40:
            dwdp[112] = 1.0*SP1;
            break;
        case 41:
            dwdp[46] = 1.0*TIMP1;
            break;
        case 42:
            dwdp[45] = 1.0*TIMP1_mRNA;
            break;
        case 43:
            dwdp[123] = 1.0*TIMP3;
            break;
        case 44:
            dwdp[122] = 1.0*TIMP3_mRNA;
            break;
        case 45:
            dwdp[67] = 1.0*cFos_P;
            break;
        case 46:
            dwdp[71] = 1.0*DUSP16*cFos_P;
            break;
        case 47:
            dwdp[13] = 1.0*cJun_P;
            break;
        case 48:
            dwdp[79] = 1.0*JAK1_P;
            break;
        case 49:
            dwdp[80] = 1.0*JAK1_P*PTPRT;
            break;
        case 50:
            dwdp[10] = 1.0*JNK_P;
            break;
        case 51:
            dwdp[11] = 1.0*DUSP16*JNK_P;
            break;
        case 52:
            dwdp[23] = 1.0*p38_P;
            break;
        case 53:
            dwdp[24] = 1.0*MKP1*p38_P;
            break;
        case 54:
            dwdp[82] = 1.0*STAT3_P_cyt;
            break;
        case 55:
            dwdp[84] = 1.0*STAT3_P_nuc;
            break;
        case 56:
            dwdp[85] = 1.0*PTPRT*STAT3_P_nuc;
            break;
        case 57:
            dwdp[83] = 1.0*PTPRT*STAT3_P_cyt;
            break;
        case 58:
            dwdp[14] = 0.5*cJun_P*(1.0*cJun_P - 1);
            break;
        case 59:
            dwdp[104] = 1.0*ADAMTS4*TIMP1;
            break;
        case 60:
            dwdp[124] = 1.0*ADAMTS4*TIMP3;
            break;
        case 61:
            dwdp[102] = 1.0*MMP13*TIMP1;
            break;
        case 62:
            dwdp[130] = 1.0*MMP13*TIMP3;
            break;
        case 63:
            dwdp[98] = 1.0*MMP1*TIMP1;
            break;
        case 64:
            dwdp[126] = 1.0*MMP1*TIMP3;
            break;
        case 65:
            dwdp[100] = 1.0*MMP3*TIMP1;
            break;
        case 66:
            dwdp[128] = 1.0*MMP3*TIMP3;
            break;
        case 67:
            dwdp[58] = 1.0*PP4*TRAF6;
            dwdp[59] = 1.0*IRAK2_TRAF6*PP4;
            break;
        case 68:
            dwdp[87] = 1.0*STAT3_nuc;
            break;
        case 69:
            dwdp[66] = 1.0*cFos*p38_P;
            break;
        case 70:
            dwdp[12] = 1.0*JNK_P*cJun;
            break;
        case 71:
            dwdp[78] = 1.0*JAK1*OSM_OSMR;
            break;
        case 72:
            dwdp[9] = 1.0*IRAK2_TRAF6*JNK;
            break;
        case 73:
            dwdp[22] = 1.0*IRAK2_TRAF6*p38;
            break;
        case 74:
            dwdp[81] = 1.0*JAK1_P*STAT3_cyt;
            break;
        case 75:
            dwdp[105] = 1.0*ADAMTS4_TIMP1;
            break;
        case 76:
            dwdp[125] = 1.0*ADAMTS4_TIMP3;
            break;
        case 77:
            dwdp[73] = 1.0*cFos_cJun;
            break;
        case 78:
            dwdp[1] = 1.0*IL1_IL1R;
            break;
        case 79:
            dwdp[3] = 1.0*IL1_IL1Ra;
            break;
        case 80:
            dwdp[6] = 1.0*IL1_IL1R_IRAK2;
            break;
        case 81:
            dwdp[99] = 1.0*MMP1_TIMP1;
            break;
        case 82:
            dwdp[103] = 1.0*MMP13_TIMP1;
            break;
        case 83:
            dwdp[131] = 1.0*MMP13_TIMP3;
            break;
        case 84:
            dwdp[127] = 1.0*MMP1_TIMP3;
            break;
        case 85:
            dwdp[101] = 1.0*MMP3_TIMP1;
            break;
        case 86:
            dwdp[129] = 1.0*MMP3_TIMP3;
            break;
        case 87:
            dwdp[75] = 1.0*OSM_OSMR;
            break;
        case 88:
            dwdp[77] = 1.0*OSM_OSMRa;
            break;
        case 89:
            dwdp[96] = 1.0*OSMR_SOCS3;
            break;
        case 90:
            dwdp[114] = 1.0*SP1_TIMP1_DNA;
            break;
        case 91:
            dwdp[8] = 1.0*IRAK2_TRAF6;
            break;
        case 92:
            dwdp[60] = 1.0*TRAF6_PP4;
            dwdp[61] = 1.0*IRAK2_TRAF6_PP4;
            break;
        case 93:
            dwdp[49] = 1.0*ADAMTS4_mRNA;
            break;
        case 94:
            dwdp[47] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 95:
            dwdp[48] = 1.0*cJun_dimer;
            break;
        case 96:
            dwdp[18] = 1.0*Source;
            break;
        case 97:
            dwdp[116] = 1.0*TIMP1_DNA;
            break;
        case 98:
            dwdp[118] = 1.0*Source;
            break;
        case 99:
            dwdp[64] = 1.0*cFos_mRNA;
            break;
        case 100:
            dwdp[62] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 101:
            dwdp[88] = 1.0*STAT3_P_nuc;
            break;
        case 102:
            dwdp[20] = 1.0*cJun_mRNA;
            break;
        case 103:
            dwdp[16] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 104:
            dwdp[17] = 1.0*cJun_dimer;
            break;
        case 105:
            dwdp[54] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 106:
            dwdp[55] = 1.0*cJun_dimer;
            break;
        case 107:
            dwdp[109] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 108:
            dwdp[68] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 109:
            dwdp[69] = 1.0*cJun_dimer;
            break;
        case 110:
            dwdp[27] = 1.0*MMP1_mRNA;
            break;
        case 111:
            dwdp[40] = 1.0*MMP13_mRNA;
            break;
        case 112:
            dwdp[38] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 113:
            dwdp[39] = 1.0*cJun_dimer;
            break;
        case 114:
            dwdp[25] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 115:
            dwdp[26] = 1.0*cJun_dimer;
            break;
        case 116:
            dwdp[34] = 1.0*MMP3_mRNA;
            break;
        case 117:
            dwdp[32] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 118:
            dwdp[33] = 1.0*cJun_dimer;
            break;
        case 119:
            dwdp[52] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 120:
            dwdp[53] = 1.0*cJun_dimer;
            break;
        case 121:
            dwdp[89] = 1.0*STAT3_P_nuc;
            break;
        case 122:
            dwdp[93] = 1.0*SOCS3_mRNA;
            break;
        case 123:
            dwdp[91] = 1.0*STAT3_P_nuc;
            break;
        case 124:
            dwdp[111] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 125:
            dwdp[44] = 1.0*TIMP1_mRNA;
            break;
        case 126:
            dwdp[117] = 1.0*TIMP1_DNA*cFos_cJun*kAP1activity;
            break;
        case 127:
            dwdp[115] = 1.0*STAT3_P_nuc*TIMP1_DNA;
            break;
        case 128:
            dwdp[121] = 1.0*TIMP3_mRNA;
            break;
        case 129:
            dwdp[119] = 1.0*cFos_cJun*kAP1activity;
            break;
        case 130:
            dwdp[120] = 1.0*STAT3_P_nuc*kAP1activity;
            break;
    }
}