#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Hui2014(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = AcanmRNA;
    y[1] = ADAMTS5;
    y[2] = AGEprod;
    y[3] = Alk1;
    y[4] = Alk1_Alk5;
    y[5] = Alk5;
    y[6] = Alk5_dimer;
    y[7] = Bax;
    y[8] = Bax_Bcl2;
    y[9] = Bax_Bcl2_Beclin;
    y[10] = Bax_Bcl2_Beclin_I;
    y[11] = Bcl2;
    y[12] = Bcl2_Beclin;
    y[13] = Bcl2_Beclin_I;
    y[14] = Beclin;
    y[15] = Beclin_I;
    y[16] = Caspase_A;
    y[17] = Caspase_I;
    y[18] = Col2mRNA;
    y[19] = DamP;
    y[20] = IkB;
    y[21] = IkB_NFkB;
    y[22] = IL1;
    y[23] = Lys_A;
    y[24] = Lys_I;
    y[25] = MMP13;
    y[26] = MMP2;
    y[27] = NatP;
    y[28] = NFkB;
    y[29] = NFkB_P;
    y[30] = p38;
    y[31] = p38_P;
    y[32] = proMMP13;
    y[33] = proMMP2;
    y[34] = RAGE;
    y[35] = ROS;
    y[36] = Runx2_A;
    y[37] = Runx2_I;
    y[38] = Smad1;
    y[39] = Smad1_P;
    y[40] = Smad1_P_Smad4;
    y[41] = Smad2;
    y[42] = Smad2_P;
    y[43] = Smad2_P_Smad4;
    y[44] = Smad4;
    y[45] = Smad7;
    y[46] = SOD;
    y[47] = Sox9;
    y[48] = Sox9_A;
    y[49] = Sox9mRNA;
    y[50] = Tgfb_A;
    y[51] = Tgfb_Alk1_Alk5;
    y[52] = Tgfb_Alk1_Alk5_Smad7;
    y[53] = Tgfb_Alk5_dimer;
    y[54] = Tgfb_Alk5_dimer_Smad7;
    y[55] = AggFrag;
    y[56] = Aggrecan;
    y[57] = Aggrecan_Collagen2;
    y[58] = ColFrag;
    y[59] = Collagen2;
    y[60] = Integrin;
    y[61] = Tgfb_I;
    y[62] = Sink;
    y[63] = Source;
    y[64] = IntegrinCount;
}