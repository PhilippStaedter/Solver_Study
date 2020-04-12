#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Proctor2013a(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ADAMTS4_mRNA;
    x_solver[1] = cFos;
    x_solver[2] = cFos_mRNA;
    x_solver[3] = cJun;
    x_solver[4] = cJun_mRNA;
    x_solver[5] = DUSP16;
    x_solver[6] = IRAK2;
    x_solver[7] = IRAK2_TRAF6;
    x_solver[8] = IRAK2_TRAF6_PP4;
    x_solver[9] = JAK1;
    x_solver[10] = JAK1_P;
    x_solver[11] = JNK;
    x_solver[12] = JNK_P;
    x_solver[13] = Matriptase;
    x_solver[14] = MKP1;
    x_solver[15] = MMP1_mRNA;
    x_solver[16] = MMP3_mRNA;
    x_solver[17] = MMP13_mRNA;
    x_solver[18] = p38;
    x_solver[19] = p38_P;
    x_solver[20] = PP4;
    x_solver[21] = proMMP1;
    x_solver[22] = proMMP3;
    x_solver[23] = proMMP13;
    x_solver[24] = PTPRT;
    x_solver[25] = SOCS3;
    x_solver[26] = SOCS3_mRNA;
    x_solver[27] = STAT3_cyt;
    x_solver[28] = STAT3_P_cyt;
    x_solver[29] = TIMP1_mRNA;
    x_solver[30] = TIMP3_mRNA;
    x_solver[31] = TRAF6;
    x_solver[32] = TRAF6_PP4;
    x_solver[33] = ADAMTS4;
    x_solver[34] = ADAMTS4_TIMP1;
    x_solver[35] = ADAMTS4_TIMP3;
    x_solver[36] = Aggrecan;
    x_solver[37] = Aggrecan_Collagen2;
    x_solver[38] = AggFrag;
    x_solver[39] = ColFrag;
    x_solver[40] = Collagen2;
    x_solver[41] = IL1;
    x_solver[42] = MMP1;
    x_solver[43] = MMP1_TIMP1;
    x_solver[44] = MMP1_TIMP3;
    x_solver[45] = MMP3;
    x_solver[46] = MMP3_TIMP1;
    x_solver[47] = MMP3_TIMP3;
    x_solver[48] = MMP13;
    x_solver[49] = MMP13_TIMP1;
    x_solver[50] = MMP13_TIMP3;
    x_solver[51] = OSM;
    x_solver[52] = TIMP1;
    x_solver[53] = TIMP3;
    x_solver[54] = IL1_IL1R;
    x_solver[55] = IL1_IL1Ra;
    x_solver[56] = IL1_IL1R_IRAK2;
    x_solver[57] = IL1R;
    x_solver[58] = IL1Ra;
    x_solver[59] = OSM_OSMR;
    x_solver[60] = OSM_OSMRa;
    x_solver[61] = OSMR_SOCS3;
    x_solver[62] = OSMR;
    x_solver[63] = OSMRa;
    x_solver[64] = cFos_cJun;
    x_solver[65] = cFos_P;
    x_solver[66] = cJun_P;
    x_solver[67] = cJun_dimer;
    x_solver[68] = SP1;
    x_solver[69] = SP1_TIMP1_DNA;
    x_solver[70] = STAT3_nuc;
    x_solver[71] = STAT3_P_nuc;
    x_solver[72] = TIMP1_DNA;
    x_solver[73] = Source;
    x_solver[74] = Sink;
}