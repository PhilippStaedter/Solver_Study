#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Hui2014(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = AcanmRNA;
    x_rdata[1] = ADAMTS5;
    x_rdata[2] = AGEprod;
    x_rdata[3] = Alk1;
    x_rdata[4] = Alk1_Alk5;
    x_rdata[5] = Alk5;
    x_rdata[6] = Alk5_dimer;
    x_rdata[7] = Bax;
    x_rdata[8] = Bax_Bcl2;
    x_rdata[9] = Bax_Bcl2_Beclin;
    x_rdata[10] = Bax_Bcl2_Beclin_I;
    x_rdata[11] = Bcl2;
    x_rdata[12] = Bcl2_Beclin;
    x_rdata[13] = Bcl2_Beclin_I;
    x_rdata[14] = Beclin;
    x_rdata[15] = Beclin_I;
    x_rdata[16] = Caspase_A;
    x_rdata[17] = Caspase_I;
    x_rdata[18] = Col2mRNA;
    x_rdata[19] = DamP;
    x_rdata[20] = IkB;
    x_rdata[21] = IkB_NFkB;
    x_rdata[22] = IL1;
    x_rdata[23] = Lys_A;
    x_rdata[24] = Lys_I;
    x_rdata[25] = MMP13;
    x_rdata[26] = MMP2;
    x_rdata[27] = NatP;
    x_rdata[28] = NFkB;
    x_rdata[29] = NFkB_P;
    x_rdata[30] = p38;
    x_rdata[31] = p38_P;
    x_rdata[32] = proMMP13;
    x_rdata[33] = proMMP2;
    x_rdata[34] = RAGE;
    x_rdata[35] = ROS;
    x_rdata[36] = Runx2_A;
    x_rdata[37] = Runx2_I;
    x_rdata[38] = Smad1;
    x_rdata[39] = Smad1_P;
    x_rdata[40] = Smad1_P_Smad4;
    x_rdata[41] = Smad2;
    x_rdata[42] = Smad2_P;
    x_rdata[43] = Smad2_P_Smad4;
    x_rdata[44] = Smad4;
    x_rdata[45] = Smad7;
    x_rdata[46] = SOD;
    x_rdata[47] = Sox9;
    x_rdata[48] = Sox9_A;
    x_rdata[49] = Sox9mRNA;
    x_rdata[50] = Tgfb_A;
    x_rdata[51] = Tgfb_Alk1_Alk5;
    x_rdata[52] = Tgfb_Alk1_Alk5_Smad7;
    x_rdata[53] = Tgfb_Alk5_dimer;
    x_rdata[54] = Tgfb_Alk5_dimer_Smad7;
    x_rdata[55] = AggFrag;
    x_rdata[56] = Aggrecan;
    x_rdata[57] = Aggrecan_Collagen2;
    x_rdata[58] = ColFrag;
    x_rdata[59] = Collagen2;
    x_rdata[60] = Integrin;
    x_rdata[61] = Tgfb_I;
    x_rdata[62] = Sink;
    x_rdata[63] = Source;
    x_rdata[64] = IntegrinCount;
}