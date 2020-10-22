#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Proctor2013a(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ADAMTS4_mRNA;
    y[1] = cFos;
    y[2] = cFos_mRNA;
    y[3] = cJun;
    y[4] = cJun_mRNA;
    y[5] = DUSP16;
    y[6] = IRAK2;
    y[7] = IRAK2_TRAF6;
    y[8] = IRAK2_TRAF6_PP4;
    y[9] = JAK1;
    y[10] = JAK1_P;
    y[11] = JNK;
    y[12] = JNK_P;
    y[13] = Matriptase;
    y[14] = MKP1;
    y[15] = MMP1_mRNA;
    y[16] = MMP3_mRNA;
    y[17] = MMP13_mRNA;
    y[18] = p38;
    y[19] = p38_P;
    y[20] = PP4;
    y[21] = proMMP1;
    y[22] = proMMP3;
    y[23] = proMMP13;
    y[24] = PTPRT;
    y[25] = SOCS3;
    y[26] = SOCS3_mRNA;
    y[27] = STAT3_cyt;
    y[28] = STAT3_P_cyt;
    y[29] = TIMP1_mRNA;
    y[30] = TIMP3_mRNA;
    y[31] = TRAF6;
    y[32] = TRAF6_PP4;
    y[33] = ADAMTS4;
    y[34] = ADAMTS4_TIMP1;
    y[35] = ADAMTS4_TIMP3;
    y[36] = Aggrecan;
    y[37] = Aggrecan_Collagen2;
    y[38] = AggFrag;
    y[39] = ColFrag;
    y[40] = Collagen2;
    y[41] = IL1;
    y[42] = MMP1;
    y[43] = MMP1_TIMP1;
    y[44] = MMP1_TIMP3;
    y[45] = MMP3;
    y[46] = MMP3_TIMP1;
    y[47] = MMP3_TIMP3;
    y[48] = MMP13;
    y[49] = MMP13_TIMP1;
    y[50] = MMP13_TIMP3;
    y[51] = OSM;
    y[52] = TIMP1;
    y[53] = TIMP3;
    y[54] = IL1_IL1R;
    y[55] = IL1_IL1Ra;
    y[56] = IL1_IL1R_IRAK2;
    y[57] = IL1R;
    y[58] = IL1Ra;
    y[59] = OSM_OSMR;
    y[60] = OSM_OSMRa;
    y[61] = OSMR_SOCS3;
    y[62] = OSMR;
    y[63] = OSMRa;
    y[64] = cFos_cJun;
    y[65] = cFos_P;
    y[66] = cJun_P;
    y[67] = cJun_dimer;
    y[68] = SP1;
    y[69] = SP1_TIMP1_DNA;
    y[70] = STAT3_nuc;
    y[71] = STAT3_P_nuc;
    y[72] = TIMP1_DNA;
    y[73] = Source;
    y[74] = Sink;
}