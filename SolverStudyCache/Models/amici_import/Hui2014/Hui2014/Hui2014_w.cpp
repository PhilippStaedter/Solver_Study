#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Hui2014(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Bax*Caspase_I*kactCasp;
    w[1] = 1.0*Beclin_I*Caspase_I*kactCaspBecI;
    w[2] = 1.0*Caspase_I*kactCaspp38*p38_P;
    w[3] = 1.0*Caspase_A*kinactCasp;
    w[4] = 1.0*Bcl2_Beclin*Caspase_A*kinactCaspBcl2;
    w[5] = 1.0*Bcl2*Caspase_A*kinactCaspBcl2;
    w[6] = 1.0*Beclin*Lys_I*kactLys;
    w[7] = 1.0*Lys_A*kinhibLys;
    w[8] = 1.0*Source*ksynBcl2;
    w[9] = 1.0*Bcl2*kdegBcl2;
    w[10] = 1.0*Bcl2*ROS*kdegBcl2ROS;
    w[11] = 1.0*Bcl2*Caspase_A*kdegBcl2Casp;
    w[12] = 1.0*Bax*Bcl2*kbinBaxBcl2;
    w[13] = 1.0*Bax_Bcl2*krelBaxBcl2;
    w[14] = 1.0*Bcl2*Beclin*kbinBcl2Beclin;
    w[15] = 1.0*Bcl2_Beclin*krelBcl2Beclin;
    w[16] = 1.0*Bcl2*Beclin_I*kbinBcl2Beclin;
    w[17] = 1.0*Bcl2_Beclin_I*krelBcl2Beclin;
    w[18] = 1.0*Beclin*kinactBec;
    w[19] = 1.0*Beclin*Caspase_A*kinactBecCasp;
    w[20] = 1.0*Bax_Bcl2*Beclin*kbinBecToBaxBcl2;
    w[21] = 1.0*Bax_Bcl2*Beclin_I*kbinBecToBaxBcl2;
    w[22] = 1.0*Bax*Bcl2_Beclin*kbinBaxToBcl2Bec;
    w[23] = 1.0*Bax*Bcl2_Beclin_I*kbinBaxToBcl2Bec;
    w[24] = 1.0*Bax_Bcl2_Beclin*krelBaxBcl2Bec;
    w[25] = 1.0*Bax_Bcl2_Beclin_I*krelBaxBcl2Bec;
    w[26] = 1.0*Bax_Bcl2_Beclin*krelBecBaxBcl2;
    w[27] = 1.0*Bax_Bcl2_Beclin_I*krelBecBaxBcl2;
    w[28] = 1.0*Source*kgenROS;
    w[29] = 1.0*ROS*kremROS;
    w[30] = 1.0*NatP*ROS*kdamNatP/(1.0*ROS + 10);
    w[31] = 1.0*DamP*Lys_A*kdegDamP;
    w[32] = 1.0*Source*kprodAGE;
    w[33] = 1.0*AGEprod*kactRAGE;
    w[34] = 1.0*RAGE*kgenROSbyRAGE;
    w[35] = 1.0*IkB_NFkB*ROS*kdegIkB;
    w[36] = 1.0*IL1*IkB_NFkB*kdegIkB;
    w[37] = 1.0*IkB*NFkB*kinactNFkB;
    w[38] = 1.0*RAGE*kinactRAGE;
    w[39] = 1.0*NFkB_P*ksynRAGE;
    w[40] = 1.0*NFkB_P*ksynIL1;
    w[41] = 1.0*IL1*kdegIL1;
    w[42] = 1.0*NFkB_P*ksynIkB;
    w[43] = 1.0*IL1*ksynMMP13;
    w[44] = 1.0*MMP13*kdegMMP13;
    w[45] = 1.0*IL1*ksynMMP2;
    w[46] = 1.0*kactMMP2*proMMP2;
    w[47] = 1.0*MMP2*kdegMMP2;
    w[48] = 1.0*IL1*ksynADAMTS5;
    w[49] = 1.0*ADAMTS5*kdegADAMTS5;
    w[50] = 1.0*ADAMTS5*Aggrecan_Collagen2*kdegAggrecan;
    w[51] = 1.0*Collagen2*MMP13*kdegCollagen;
    w[52] = 1.0*DamP*kgenROSbyDamP;
    w[53] = 1.0*Source*ksynNatP;
    w[54] = 1.0*NFkB_P*ksynSOD;
    w[55] = 1.0*SOD*kdegSOD;
    w[56] = 1.0*ROS*SOD*kremROSbySOD;
    w[57] = 1.0*IL1*kphosp38*p38;
    w[58] = 1.0*ROS*kphosp38ROS*p38;
    w[59] = 1.0*kdephosp38*p38_P;
    w[60] = 1.0*NFkB*kphosNFkB*p38_P;
    w[61] = 1.0*NFkB_P*kdephosNFkB;
    w[62] = 1.0*kgenROSbyp38*p38_P;
    w[63] = 1.0*Lys_A*ROS*kdamLys/(1.0*ROS + 10);
    w[64] = 1.0*Source*kactIntegrin;
    w[65] = 1.0*Integrin*kinactIntegrin;
    w[66] = 1.0*Source*ksynAlk5;
    w[67] = 1.0*Integrin*Tgfb_I*kactTgfbIntegrin;
    w[68] = 1.0*MMP2*Tgfb_I*kactTgfbMMP2;
    w[69] = 1.0*Tgfb_A*kinactTgfb;
    w[70] = 0.5*Alk5*kdimerAlk5*(1.0*Alk5 - 1);
    w[71] = 1.0*Alk5_dimer*kdedimerAlk5;
    w[72] = 1.0*Alk1*Alk5*kbinAlk1Alk5;
    w[73] = 1.0*Alk1_Alk5*krelAlk1Alk5;
    w[74] = 1.0*Alk5_dimer*Tgfb_A*kbinTgfbAlk5;
    w[75] = 1.0*Tgfb_Alk5_dimer*krelTgfbAlk5;
    w[76] = 1.0*Smad7*Tgfb_Alk5_dimer*kbinSmad7Alk5;
    w[77] = 1.0*Tgfb_Alk5_dimer_Smad7*krelSmad7Alk5;
    w[78] = 1.0*Tgfb_Alk5_dimer_Smad7*kdegSmad7Alk5;
    w[79] = 1.0*Alk1_Alk5*Tgfb_A*kbinTgfbAlk1;
    w[80] = 1.0*Tgfb_Alk1_Alk5*krelTgfbAlk1;
    w[81] = 1.0*Smad2*Tgfb_Alk5_dimer*kphosSmad2;
    w[82] = 1.0*Smad2_P*Smad4*kbinSmad2Smad4;
    w[83] = 1.0*Smad2_P_Smad4*krelSmad2Smad4;
    w[84] = 1.0*Smad2_P*kdephosSmad2;
    w[85] = 1.0*Smad2_P_Smad4*ksynSmad7;
    w[86] = 1.0*Smad2_P_Smad4*Sox9*kactSox9;
    w[87] = 1.0*Sox9_A*kinactSox9;
    w[88] = 1.0*Source*ksynSox9mRNA;
    w[89] = 1.0*Sox9_A*ksynSox9mRNASox9A;
    w[90] = 1.0*Sox9mRNA*kdegSox9mRNA;
    w[91] = 1.0*Sox9mRNA*ksynSox9;
    w[92] = 1.0*Sox9*kdegSox9;
    w[93] = 1.0*Sox9_A*ksynCol2mRNASox9A;
    w[94] = 1.0*Smad2_P_Smad4*ksynCol2mRNASmad;
    w[95] = 1.0*Col2mRNA*kdegCol2mRNA;
    w[96] = 1.0*Col2mRNA*ksynCol2;
    w[97] = 1.0*Sox9_A*ksynAcanmRNASox9A;
    w[98] = 1.0*AcanmRNA*kdegAcanmRNA;
    w[99] = 1.0*AcanmRNA*ksynAggrecan;
    w[100] = 1.0*Aggrecan*Collagen2*kbinAggrecanCollagen2;
    w[101] = 1.0*Runx2_A*Smad2_P_Smad4*kinactRunx2;
    w[102] = 1.0*Alk5*kdegAlk5;
    w[103] = 1.0*Smad1*Tgfb_Alk1_Alk5*kphosSmad1;
    w[104] = 1.0*Smad1_P*kdephosSmad1;
    w[105] = 1.0*Smad1_P*Smad7*kdephosSmad1Smad7;
    w[106] = 1.0*Smad1_P*Smad4*kbinSmad1Smad4;
    w[107] = 1.0*Smad1_P_Smad4*krelSmad1Smad4;
    w[108] = 1.0*Runx2_I*Smad1_P_Smad4*kactRunx2;
    w[109] = 1.0*Runx2_A*ksynMMP13Runx2;
    w[110] = 1.0*Source*ksynAlk1;
    w[111] = 1.0*Alk1*kdegAlk1;
    w[112] = 1.0*kactMMP13*proMMP13;
    w[113] = 1.0*Smad7*Tgfb_Alk1_Alk5*kbinSmad7Alk1;
    w[114] = 1.0*Tgfb_Alk1_Alk5_Smad7*krelSmad7Alk1;
    w[115] = 1.0*Tgfb_Alk1_Alk5_Smad7*kdegSmad7Alk1;
    w[116] = 1.0*Smad7*kdegSmad7;
}