#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Proctor2013a(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*ksynADAMTS4;
    dwdx[1] = 1.0*kdegADAMTS4mRNA;
    dwdx[2] = 1.0*kdegcFos;
    dwdx[3] = 1.0*kphoscFos*p38_P;
    dwdx[4] = 1.0*kdegcFosmRNA;
    dwdx[5] = 1.0*ksyncFos;
    dwdx[6] = 1.0*JNK_P*kphoscJun;
    dwdx[7] = 1.0*kdegcJun;
    dwdx[8] = 1.0*kdegcJunmRNA;
    dwdx[9] = 1.0*ksyncJun;
    dwdx[10] = 1.0*JNK_P*kdephosJNKDUSP16;
    dwdx[11] = 1.0*kdegDUSP16;
    dwdx[12] = 1.0*cFos_P*kdephoscFosDUSP16;
    dwdx[13] = 1.0*IL1_IL1R*kbinIRAK2;
    dwdx[14] = 1.0*krelTRAF6;
    dwdx[15] = 1.0*JNK*kphosJNK;
    dwdx[16] = 1.0*kphosp38*p38;
    dwdx[17] = 1.0*PP4*kinhibTRAF6;
    dwdx[18] = 1.0*krelTRAF6PP4;
    dwdx[19] = 1.0*OSM_OSMR*kphosJAK1;
    dwdx[20] = 1.0*kdephosJAK1;
    dwdx[21] = 1.0*PTPRT*kdephosJAK1PTPRT;
    dwdx[22] = 1.0*STAT3_cyt*kphosSTAT3;
    dwdx[23] = 1.0*IRAK2_TRAF6*kphosJNK;
    dwdx[24] = 1.0*kdephosJNK;
    dwdx[25] = 1.0*DUSP16*kdephosJNKDUSP16;
    dwdx[26] = 1.0*cJun*kphoscJun;
    dwdx[27] = 1.0*kactMMP1mat*proMMP1;
    dwdx[28] = 1.0*kactMMP3mat*proMMP3;
    dwdx[29] = 1.0*kdegMatriptase;
    dwdx[30] = 1.0*kdephosp38MKP1*p38_P;
    dwdx[31] = 1.0*kdegMKP1;
    dwdx[32] = 1.0*ksynMMP1;
    dwdx[33] = 1.0*kdegMMP1mRNA;
    dwdx[34] = 1.0*ksynMMP3;
    dwdx[35] = 1.0*kdegMMP3mRNA;
    dwdx[36] = 1.0*ksynMMP13;
    dwdx[37] = 1.0*kdegMMP13mRNA;
    dwdx[38] = 1.0*IRAK2_TRAF6*kphosp38;
    dwdx[39] = 1.0*kdephosp38;
    dwdx[40] = 1.0*MKP1*kdephosp38MKP1;
    dwdx[41] = 1.0*cFos*kphoscFos;
    dwdx[42] = 1.0*kdegPP4;
    dwdx[43] = 1.0*TRAF6*kinhibTRAF6;
    dwdx[44] = 1.0*IRAK2_TRAF6*kinhibTRAF6;
    dwdx[45] = 1.0*Matriptase*kactMMP1mat;
    dwdx[46] = 1.0*MMP3*kactMMP1mmp3;
    dwdx[47] = 1.0*Matriptase*kactMMP3mat;
    dwdx[48] = 1.0*MMP3*kactMMP13mmp3;
    dwdx[49] = 1.0*JAK1_P*kdephosJAK1PTPRT;
    dwdx[50] = 1.0*STAT3_P_cyt*kdephosSTAT3PTPRT;
    dwdx[51] = 1.0*STAT3_P_nuc*kdephosSTAT3nucPTPRT;
    dwdx[52] = 1.0*kdegPTPRT;
    dwdx[53] = 1.0*kdegSOCS3;
    dwdx[54] = 1.0*OSMR*kbinSOCS3OSMR;
    dwdx[55] = 1.0*kdegSOCS3mRNA;
    dwdx[56] = 1.0*ksynSOCS3;
    dwdx[57] = 1.0*JAK1_P*kphosSTAT3;
    dwdx[58] = 1.0*kdephosSTAT3;
    dwdx[59] = 1.0*PTPRT*kdephosSTAT3PTPRT;
    dwdx[60] = 1.0*kcyt2nucSTAT3;
    dwdx[61] = 1.0*ksynTIMP1;
    dwdx[62] = 1.0*kdegTIMP1mRNA;
    dwdx[63] = 1.0*ksynTIMP3;
    dwdx[64] = 1.0*kdegTIMP3mRNA;
    dwdx[65] = 1.0*IL1_IL1R_IRAK2*kbinTRAF6;
    dwdx[66] = 1.0*PP4*kinhibTRAF6;
    dwdx[67] = 1.0*krelTRAF6PP4;
    dwdx[68] = 1.0*kdegADAMTS4;
    dwdx[69] = 1.0*TIMP1*kinhibADAMTS4TIMP1;
    dwdx[70] = 1.0*Aggrecan_Collagen2*kdegAggrecan;
    dwdx[71] = 1.0*TIMP3*kinhibADAMTS4TIMP3;
    dwdx[72] = 1.0*krelADAMTS4TIMP1;
    dwdx[73] = 1.0*krelADAMTS4TIMP3;
    dwdx[74] = 1.0*ADAMTS4*kdegAggrecan;
    dwdx[75] = 1.0*MMP1*kdegCollagen2mmp1;
    dwdx[76] = 1.0*MMP13*kdegCollagen2mmp13;
    dwdx[77] = 1.0*IL1R*kbinIL1IL1R;
    dwdx[78] = 1.0*IL1Ra*kbinIL1IL1Ra;
    dwdx[79] = 1.0*kdegIL1;
    dwdx[80] = 1.0*kdegMMP1;
    dwdx[81] = 1.0*TIMP1*kinhibMMP1TIMP1;
    dwdx[82] = 1.0*Collagen2*kdegCollagen2mmp1;
    dwdx[83] = 1.0*TIMP3*kinhibMMP1TIMP3;
    dwdx[84] = 1.0*krelMMP1;
    dwdx[85] = 1.0*krelMMP1TIMP3;
    dwdx[86] = 1.0*kactMMP1mmp3*proMMP1;
    dwdx[87] = 1.0*kdegMMP3;
    dwdx[88] = 1.0*kactMMP13mmp3*proMMP13;
    dwdx[89] = 1.0*TIMP1*kinhibMMP3TIMP1;
    dwdx[90] = 1.0*TIMP3*kinhibMMP3TIMP3;
    dwdx[91] = 1.0*krelMMP3;
    dwdx[92] = 1.0*krelMMP3TIMP3;
    dwdx[93] = 1.0*kdegMMP13;
    dwdx[94] = 1.0*TIMP1*kinhibMMP13TIMP1;
    dwdx[95] = 1.0*Collagen2*kdegCollagen2mmp13;
    dwdx[96] = 1.0*TIMP3*kinhibMMP13TIMP3;
    dwdx[97] = 1.0*krelMMP13;
    dwdx[98] = 1.0*krelMMP13TIMP3;
    dwdx[99] = 1.0*OSMR*kbinOSMOSMR;
    dwdx[100] = 1.0*OSMRa*kbinOSMOSMRa;
    dwdx[101] = 1.0*kdegOSM;
    dwdx[102] = 1.0*kdegTIMP1;
    dwdx[103] = 1.0*MMP1*kinhibMMP1TIMP1;
    dwdx[104] = 1.0*MMP3*kinhibMMP3TIMP1;
    dwdx[105] = 1.0*MMP13*kinhibMMP13TIMP1;
    dwdx[106] = 1.0*ADAMTS4*kinhibADAMTS4TIMP1;
    dwdx[107] = 1.0*kdegTIMP3;
    dwdx[108] = 1.0*ADAMTS4*kinhibADAMTS4TIMP3;
    dwdx[109] = 1.0*MMP1*kinhibMMP1TIMP3;
    dwdx[110] = 1.0*MMP3*kinhibMMP3TIMP3;
    dwdx[111] = 1.0*MMP13*kinhibMMP13TIMP3;
    dwdx[112] = 1.0*krelIL1IL1R;
    dwdx[113] = 1.0*IRAK2*kbinIRAK2;
    dwdx[114] = 1.0*krelIL1IL1Ra;
    dwdx[115] = 1.0*krelIRAK2;
    dwdx[116] = 1.0*TRAF6*kbinTRAF6;
    dwdx[117] = 1.0*IL1*kbinIL1IL1R;
    dwdx[118] = 1.0*IL1*kbinIL1IL1Ra;
    dwdx[119] = 1.0*krelOSMOSMR;
    dwdx[120] = 1.0*JAK1*kphosJAK1;
    dwdx[121] = 1.0*krelOSMOSMRa;
    dwdx[122] = 1.0*krelSOCS3OSMR;
    dwdx[123] = 1.0*OSM*kbinOSMOSMR;
    dwdx[124] = 1.0*SOCS3*kbinSOCS3OSMR;
    dwdx[125] = 1.0*OSM*kbinOSMOSMRa;
    dwdx[126] = 1.0*kAP1activity*ksyncJunmRNA;
    dwdx[127] = 1.0*kAP1activity*ksynMMP1mRNA;
    dwdx[128] = 1.0*kAP1activity*ksynMMP3mRNA;
    dwdx[129] = 1.0*kAP1activity*ksynMMP13mRNA;
    dwdx[130] = 1.0*kAP1activity*ksynADAMTS4mRNA;
    dwdx[131] = 1.0*kAP1activity*ksynPP4;
    dwdx[132] = 1.0*kAP1activity*ksynDUSP16;
    dwdx[133] = 1.0*kAP1activity*ksyncFosmRNA;
    dwdx[134] = 1.0*kAP1activity*ksynMKP1;
    dwdx[135] = 1.0*krelcFoscJun;
    dwdx[136] = 1.0*kAP1activity*ksynMatriptase;
    dwdx[137] = 1.0*kAP1activity*ksynSP1;
    dwdx[138] = 1.0*TIMP1_DNA*kAP1activity*ksynTIMP1mRNA;
    dwdx[139] = 1.0*kAP1activity*ksynTIMP3mRNA;
    dwdx[140] = 1.0*kdephoscFos;
    dwdx[141] = 1.0*DUSP16*kdephoscFosDUSP16;
    dwdx[142] = 1.0*cJun_P*kbincFoscJun;
    dwdx[143] = 1.0*kdephoscJun;
    dwdx[144] = 0.5*cJun_P*kdimercJun + 0.5*kdimercJun*(1.0*cJun_P - 1);
    dwdx[145] = 1.0*cFos_P*kbincFoscJun;
    dwdx[146] = 1.0*kdedimercJun;
    dwdx[147] = 1.0*ksyncJunmRNAcJun;
    dwdx[148] = 1.0*ksynMMP1mRNAcJun;
    dwdx[149] = 1.0*ksynMMP3mRNAcJun;
    dwdx[150] = 1.0*ksynMMP13mRNAcJun;
    dwdx[151] = 1.0*ksynADAMTS4mRNAcJun;
    dwdx[152] = 1.0*ksynPP4cJun;
    dwdx[153] = 1.0*ksynDUSP16cJun;
    dwdx[154] = 1.0*ksynMKP1cJun;
    dwdx[155] = 1.0*kdegSP1;
    dwdx[156] = 1.0*TIMP1_DNA*kbinSP1TIMP1DNA;
    dwdx[157] = 1.0*krelSP1TIMP1DNA;
    dwdx[158] = 1.0*knuc2cytSTAT3;
    dwdx[159] = 1.0*kdephosSTAT3nuc;
    dwdx[160] = 1.0*PTPRT*kdephosSTAT3nucPTPRT;
    dwdx[161] = 1.0*ksyncFosmRNAStat3;
    dwdx[162] = 1.0*ksynPTPRT;
    dwdx[163] = 1.0*ksynSOCS3mRNA;
    dwdx[164] = 1.0*TIMP1_DNA*ksynTIMP1mRNAStat3;
    dwdx[165] = 1.0*kAP1activity*ksynTIMP3mRNAStat3;
    dwdx[166] = 1.0*SP1*kbinSP1TIMP1DNA;
    dwdx[167] = 1.0*STAT3_P_nuc*ksynTIMP1mRNAStat3;
    dwdx[168] = 1.0*ksynbasalTIMP1mRNA;
    dwdx[169] = 1.0*cFos_cJun*kAP1activity*ksynTIMP1mRNA;
    dwdx[170] = 1.0*ksynbasalcJunmRNA;
    dwdx[171] = 1.0*ksynbasalTIMP3mRNA;
}