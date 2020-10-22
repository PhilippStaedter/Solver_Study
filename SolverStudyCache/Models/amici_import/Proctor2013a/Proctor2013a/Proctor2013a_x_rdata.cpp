#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Proctor2013a(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ADAMTS4_mRNA;
    x_rdata[1] = cFos;
    x_rdata[2] = cFos_mRNA;
    x_rdata[3] = cJun;
    x_rdata[4] = cJun_mRNA;
    x_rdata[5] = DUSP16;
    x_rdata[6] = IRAK2;
    x_rdata[7] = IRAK2_TRAF6;
    x_rdata[8] = IRAK2_TRAF6_PP4;
    x_rdata[9] = JAK1;
    x_rdata[10] = JAK1_P;
    x_rdata[11] = JNK;
    x_rdata[12] = JNK_P;
    x_rdata[13] = Matriptase;
    x_rdata[14] = MKP1;
    x_rdata[15] = MMP1_mRNA;
    x_rdata[16] = MMP3_mRNA;
    x_rdata[17] = MMP13_mRNA;
    x_rdata[18] = p38;
    x_rdata[19] = p38_P;
    x_rdata[20] = PP4;
    x_rdata[21] = proMMP1;
    x_rdata[22] = proMMP3;
    x_rdata[23] = proMMP13;
    x_rdata[24] = PTPRT;
    x_rdata[25] = SOCS3;
    x_rdata[26] = SOCS3_mRNA;
    x_rdata[27] = STAT3_cyt;
    x_rdata[28] = STAT3_P_cyt;
    x_rdata[29] = TIMP1_mRNA;
    x_rdata[30] = TIMP3_mRNA;
    x_rdata[31] = TRAF6;
    x_rdata[32] = TRAF6_PP4;
    x_rdata[33] = ADAMTS4;
    x_rdata[34] = ADAMTS4_TIMP1;
    x_rdata[35] = ADAMTS4_TIMP3;
    x_rdata[36] = Aggrecan;
    x_rdata[37] = Aggrecan_Collagen2;
    x_rdata[38] = AggFrag;
    x_rdata[39] = ColFrag;
    x_rdata[40] = Collagen2;
    x_rdata[41] = IL1;
    x_rdata[42] = MMP1;
    x_rdata[43] = MMP1_TIMP1;
    x_rdata[44] = MMP1_TIMP3;
    x_rdata[45] = MMP3;
    x_rdata[46] = MMP3_TIMP1;
    x_rdata[47] = MMP3_TIMP3;
    x_rdata[48] = MMP13;
    x_rdata[49] = MMP13_TIMP1;
    x_rdata[50] = MMP13_TIMP3;
    x_rdata[51] = OSM;
    x_rdata[52] = TIMP1;
    x_rdata[53] = TIMP3;
    x_rdata[54] = IL1_IL1R;
    x_rdata[55] = IL1_IL1Ra;
    x_rdata[56] = IL1_IL1R_IRAK2;
    x_rdata[57] = IL1R;
    x_rdata[58] = IL1Ra;
    x_rdata[59] = OSM_OSMR;
    x_rdata[60] = OSM_OSMRa;
    x_rdata[61] = OSMR_SOCS3;
    x_rdata[62] = OSMR;
    x_rdata[63] = OSMRa;
    x_rdata[64] = cFos_cJun;
    x_rdata[65] = cFos_P;
    x_rdata[66] = cJun_P;
    x_rdata[67] = cJun_dimer;
    x_rdata[68] = SP1;
    x_rdata[69] = SP1_TIMP1_DNA;
    x_rdata[70] = STAT3_nuc;
    x_rdata[71] = STAT3_P_nuc;
    x_rdata[72] = TIMP1_DNA;
    x_rdata[73] = Source;
    x_rdata[74] = Sink;
}