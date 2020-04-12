#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Hui2014(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = AcanmRNA;
    x_solver[1] = ADAMTS5;
    x_solver[2] = AGEprod;
    x_solver[3] = Alk1;
    x_solver[4] = Alk1_Alk5;
    x_solver[5] = Alk5;
    x_solver[6] = Alk5_dimer;
    x_solver[7] = Bax;
    x_solver[8] = Bax_Bcl2;
    x_solver[9] = Bax_Bcl2_Beclin;
    x_solver[10] = Bax_Bcl2_Beclin_I;
    x_solver[11] = Bcl2;
    x_solver[12] = Bcl2_Beclin;
    x_solver[13] = Bcl2_Beclin_I;
    x_solver[14] = Beclin;
    x_solver[15] = Beclin_I;
    x_solver[16] = Caspase_A;
    x_solver[17] = Caspase_I;
    x_solver[18] = Col2mRNA;
    x_solver[19] = DamP;
    x_solver[20] = IkB;
    x_solver[21] = IkB_NFkB;
    x_solver[22] = IL1;
    x_solver[23] = Lys_A;
    x_solver[24] = Lys_I;
    x_solver[25] = MMP13;
    x_solver[26] = MMP2;
    x_solver[27] = NatP;
    x_solver[28] = NFkB;
    x_solver[29] = NFkB_P;
    x_solver[30] = p38;
    x_solver[31] = p38_P;
    x_solver[32] = proMMP13;
    x_solver[33] = proMMP2;
    x_solver[34] = RAGE;
    x_solver[35] = ROS;
    x_solver[36] = Runx2_A;
    x_solver[37] = Runx2_I;
    x_solver[38] = Smad1;
    x_solver[39] = Smad1_P;
    x_solver[40] = Smad1_P_Smad4;
    x_solver[41] = Smad2;
    x_solver[42] = Smad2_P;
    x_solver[43] = Smad2_P_Smad4;
    x_solver[44] = Smad4;
    x_solver[45] = Smad7;
    x_solver[46] = SOD;
    x_solver[47] = Sox9;
    x_solver[48] = Sox9_A;
    x_solver[49] = Sox9mRNA;
    x_solver[50] = Tgfb_A;
    x_solver[51] = Tgfb_Alk1_Alk5;
    x_solver[52] = Tgfb_Alk1_Alk5_Smad7;
    x_solver[53] = Tgfb_Alk5_dimer;
    x_solver[54] = Tgfb_Alk5_dimer_Smad7;
    x_solver[55] = AggFrag;
    x_solver[56] = Aggrecan;
    x_solver[57] = Aggrecan_Collagen2;
    x_solver[58] = ColFrag;
    x_solver[59] = Collagen2;
    x_solver[60] = Integrin;
    x_solver[61] = Tgfb_I;
    x_solver[62] = Sink;
    x_solver[63] = Source;
    x_solver[64] = IntegrinCount;
}