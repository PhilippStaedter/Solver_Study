#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Hui2014(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*Bax*Caspase_I;
            break;
        case 1:
            dwdp[1] = 1.0*Beclin_I*Caspase_I;
            break;
        case 2:
            dwdp[2] = 1.0*Caspase_I*p38_P;
            break;
        case 3:
            dwdp[64] = 1.0*Source;
            break;
        case 4:
            dwdp[6] = 1.0*Beclin*Lys_I;
            break;
        case 5:
            dwdp[112] = 1.0*proMMP13;
            break;
        case 6:
            dwdp[46] = 1.0*proMMP2;
            break;
        case 8:
            dwdp[33] = 1.0*AGEprod;
            break;
        case 9:
            dwdp[108] = 1.0*Runx2_I*Smad1_P_Smad4;
            break;
        case 10:
            dwdp[86] = 1.0*Smad2_P_Smad4*Sox9;
            break;
        case 11:
            dwdp[67] = 1.0*Integrin*Tgfb_I;
            break;
        case 12:
            dwdp[68] = 1.0*MMP2*Tgfb_I;
            break;
        case 13:
            dwdp[100] = 1.0*Aggrecan*Collagen2;
            break;
        case 14:
            dwdp[72] = 1.0*Alk1*Alk5;
            break;
        case 15:
            dwdp[12] = 1.0*Bax*Bcl2;
            break;
        case 16:
            dwdp[22] = 1.0*Bax*Bcl2_Beclin;
            dwdp[23] = 1.0*Bax*Bcl2_Beclin_I;
            break;
        case 17:
            dwdp[14] = 1.0*Bcl2*Beclin;
            dwdp[16] = 1.0*Bcl2*Beclin_I;
            break;
        case 18:
            dwdp[20] = 1.0*Bax_Bcl2*Beclin;
            dwdp[21] = 1.0*Bax_Bcl2*Beclin_I;
            break;
        case 19:
            dwdp[106] = 1.0*Smad1_P*Smad4;
            break;
        case 20:
            dwdp[82] = 1.0*Smad2_P*Smad4;
            break;
        case 21:
            dwdp[113] = 1.0*Smad7*Tgfb_Alk1_Alk5;
            break;
        case 22:
            dwdp[76] = 1.0*Smad7*Tgfb_Alk5_dimer;
            break;
        case 23:
            dwdp[79] = 1.0*Alk1_Alk5*Tgfb_A;
            break;
        case 24:
            dwdp[74] = 1.0*Alk5_dimer*Tgfb_A;
            break;
        case 25:
            dwdp[63] = 1.0*Lys_A*ROS/(1.0*ROS + 10);
            break;
        case 26:
            dwdp[30] = 1.0*NatP*ROS/(1.0*ROS + 10);
            break;
        case 27:
            dwdp[71] = 1.0*Alk5_dimer;
            break;
        case 28:
            dwdp[98] = 1.0*AcanmRNA;
            break;
        case 29:
            dwdp[49] = 1.0*ADAMTS5;
            break;
        case 30:
            dwdp[50] = 1.0*ADAMTS5*Aggrecan_Collagen2;
            break;
        case 31:
            dwdp[111] = 1.0*Alk1;
            break;
        case 32:
            dwdp[102] = 1.0*Alk5;
            break;
        case 33:
            dwdp[9] = 1.0*Bcl2;
            break;
        case 34:
            dwdp[11] = 1.0*Bcl2*Caspase_A;
            break;
        case 35:
            dwdp[10] = 1.0*Bcl2*ROS;
            break;
        case 36:
            dwdp[95] = 1.0*Col2mRNA;
            break;
        case 37:
            dwdp[51] = 1.0*Collagen2*MMP13;
            break;
        case 38:
            dwdp[31] = 1.0*DamP*Lys_A;
            break;
        case 39:
            dwdp[35] = 1.0*IkB_NFkB*ROS;
            dwdp[36] = 1.0*IL1*IkB_NFkB;
            break;
        case 40:
            dwdp[41] = 1.0*IL1;
            break;
        case 41:
            dwdp[44] = 1.0*MMP13;
            break;
        case 42:
            dwdp[47] = 1.0*MMP2;
            break;
        case 43:
            dwdp[116] = 1.0*Smad7;
            break;
        case 44:
            dwdp[115] = 1.0*Tgfb_Alk1_Alk5_Smad7;
            break;
        case 45:
            dwdp[78] = 1.0*Tgfb_Alk5_dimer_Smad7;
            break;
        case 46:
            dwdp[55] = 1.0*SOD;
            break;
        case 47:
            dwdp[92] = 1.0*Sox9;
            break;
        case 48:
            dwdp[90] = 1.0*Sox9mRNA;
            break;
        case 50:
            dwdp[61] = 1.0*NFkB_P;
            break;
        case 51:
            dwdp[59] = 1.0*p38_P;
            break;
        case 52:
            dwdp[104] = 1.0*Smad1_P;
            break;
        case 53:
            dwdp[105] = 1.0*Smad1_P*Smad7;
            break;
        case 54:
            dwdp[84] = 1.0*Smad2_P;
            break;
        case 55:
            dwdp[70] = 0.5*Alk5*(1.0*Alk5 - 1);
            break;
        case 56:
            dwdp[28] = 1.0*Source;
            break;
        case 57:
            dwdp[52] = 1.0*DamP;
            break;
        case 58:
            dwdp[62] = 1.0*p38_P;
            break;
        case 59:
            dwdp[34] = 1.0*RAGE;
            break;
        case 60:
            dwdp[18] = 1.0*Beclin;
            break;
        case 61:
            dwdp[19] = 1.0*Beclin*Caspase_A;
            break;
        case 62:
            dwdp[3] = 1.0*Caspase_A;
            break;
        case 63:
            dwdp[4] = 1.0*Bcl2_Beclin*Caspase_A;
            dwdp[5] = 1.0*Bcl2*Caspase_A;
            break;
        case 64:
            dwdp[65] = 1.0*Integrin;
            break;
        case 65:
            dwdp[37] = 1.0*IkB*NFkB;
            break;
        case 66:
            dwdp[38] = 1.0*RAGE;
            break;
        case 67:
            dwdp[101] = 1.0*Runx2_A*Smad2_P_Smad4;
            break;
        case 68:
            dwdp[87] = 1.0*Sox9_A;
            break;
        case 69:
            dwdp[69] = 1.0*Tgfb_A;
            break;
        case 70:
            dwdp[7] = 1.0*Lys_A;
            break;
        case 71:
            dwdp[60] = 1.0*NFkB*p38_P;
            break;
        case 72:
            dwdp[57] = 1.0*IL1*p38;
            break;
        case 73:
            dwdp[58] = 1.0*ROS*p38;
            break;
        case 74:
            dwdp[103] = 1.0*Smad1*Tgfb_Alk1_Alk5;
            break;
        case 75:
            dwdp[81] = 1.0*Smad2*Tgfb_Alk5_dimer;
            break;
        case 76:
            dwdp[32] = 1.0*Source;
            break;
        case 77:
            dwdp[73] = 1.0*Alk1_Alk5;
            break;
        case 78:
            dwdp[13] = 1.0*Bax_Bcl2;
            break;
        case 79:
            dwdp[24] = 1.0*Bax_Bcl2_Beclin;
            dwdp[25] = 1.0*Bax_Bcl2_Beclin_I;
            break;
        case 80:
            dwdp[15] = 1.0*Bcl2_Beclin;
            dwdp[17] = 1.0*Bcl2_Beclin_I;
            break;
        case 81:
            dwdp[26] = 1.0*Bax_Bcl2_Beclin;
            dwdp[27] = 1.0*Bax_Bcl2_Beclin_I;
            break;
        case 82:
            dwdp[107] = 1.0*Smad1_P_Smad4;
            break;
        case 83:
            dwdp[83] = 1.0*Smad2_P_Smad4;
            break;
        case 84:
            dwdp[114] = 1.0*Tgfb_Alk1_Alk5_Smad7;
            break;
        case 85:
            dwdp[77] = 1.0*Tgfb_Alk5_dimer_Smad7;
            break;
        case 86:
            dwdp[80] = 1.0*Tgfb_Alk1_Alk5;
            break;
        case 87:
            dwdp[75] = 1.0*Tgfb_Alk5_dimer;
            break;
        case 88:
            dwdp[29] = 1.0*ROS;
            break;
        case 89:
            dwdp[56] = 1.0*ROS*SOD;
            break;
        case 91:
            dwdp[97] = 1.0*Sox9_A;
            break;
        case 92:
            dwdp[48] = 1.0*IL1;
            break;
        case 93:
            dwdp[99] = 1.0*AcanmRNA;
            break;
        case 94:
            dwdp[110] = 1.0*Source;
            break;
        case 95:
            dwdp[66] = 1.0*Source;
            break;
        case 96:
            dwdp[8] = 1.0*Source;
            break;
        case 97:
            dwdp[96] = 1.0*Col2mRNA;
            break;
        case 99:
            dwdp[94] = 1.0*Smad2_P_Smad4;
            break;
        case 100:
            dwdp[93] = 1.0*Sox9_A;
            break;
        case 101:
            dwdp[42] = 1.0*NFkB_P;
            break;
        case 102:
            dwdp[40] = 1.0*NFkB_P;
            break;
        case 103:
            dwdp[43] = 1.0*IL1;
            break;
        case 104:
            dwdp[109] = 1.0*Runx2_A;
            break;
        case 105:
            dwdp[45] = 1.0*IL1;
            break;
        case 106:
            dwdp[53] = 1.0*Source;
            break;
        case 107:
            dwdp[39] = 1.0*NFkB_P;
            break;
        case 108:
            dwdp[85] = 1.0*Smad2_P_Smad4;
            break;
        case 109:
            dwdp[54] = 1.0*NFkB_P;
            break;
        case 110:
            dwdp[91] = 1.0*Sox9mRNA;
            break;
        case 111:
            dwdp[88] = 1.0*Source;
            break;
        case 112:
            dwdp[89] = 1.0*Sox9_A;
            break;
    }
}