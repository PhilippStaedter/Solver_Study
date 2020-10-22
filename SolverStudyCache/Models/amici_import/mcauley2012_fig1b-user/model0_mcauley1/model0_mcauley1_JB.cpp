#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_mcauley1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[11] = -1.0*dwdx0;
    JB[35] = 1.0*dwdx1;
    JB[36] = 1.0*dwdx1;
    JB[57] = -1.0*dwdx1;
    JB[69] = -1.0*dwdx2 + 1.0*dwdx6 - 1.0*dwdx8;
    JB[70] = 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7;
    JB[85] = -1.0*dwdx3;
    JB[88] = -1.0*dwdx4;
    JB[89] = -1.0*dwdx5;
    JB[91] = -1.0*dwdx6;
    JB[140] = 1.0*dwdx9;
    JB[167] = -1.0*dwdx9;
    JB[172] = 1.0*dwdx11;
    JB[174] = -1.0*dwdx10;
    JB[190] = -1.0*dwdx11;
    JB[201] = 1.0*dwdx10;
    JB[206] = -1.0*dwdx13;
    JB[208] = 1.0*dwdx12;
    JB[224] = 1.0*dwdx13;
    JB[235] = -1.0*dwdx12;
    JB[280] = 1.0*dwdx14 + 1.0*dwdx15;
    JB[285] = -1.0*dwdx15;
    JB[303] = -1.0*dwdx14;
    JB[315] = 1.0*dwdx16;
    JB[318] = -1.0*dwdx16;
    JB[321] = 1.0*dwdx17;
    JB[337] = -1.0*dwdx17;
    JB[349] = -1.0*dwdx18;
    JB[385] = -1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21;
    JB[405] = -1.0*dwdx20;
    JB[406] = -1.0*dwdx21;
    JB[455] = 1.0*dwdx22 + 1.0*dwdx23;
    JB[457] = -1.0*dwdx23;
    JB[473] = -1.0*dwdx22;
    JB[484] = 1.0*dwdx24;
    JB[489] = -1.0*dwdx24;
    JB[512] = -1.0*dwdx27 - 1.0*dwdx28;
    JB[525] = 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27 + 1.0*dwdx28;
    JB[541] = -1.0*dwdx25 - 1.0*dwdx26;
    JB[557] = 1.0*dwdx29;
    JB[559] = -1.0*dwdx29;
    JB[580] = -1.0*dwdx30;
    JB[593] = 1.0*dwdx30;
    JB[595] = 1.0*dwdx31;
    JB[597] = -1.0*dwdx31;
    JB[629] = -1.0*dwdx32;
    JB[682] = -1.0*dwdx33;
    JB[700] = 1.0*dwdx33;
    JB[790] = -1.0*dwdx34;
    JB[797] = -1.0*dwdx35;
    JB[805] = 1.0*dwdx34 + 1.0*dwdx35 + 1.0*dwdx36;
    JB[813] = -1.0*dwdx36;
    JB[817] = 1.0*dwdx37;
    JB[818] = 1.0*dwdx37;
    JB[839] = -1.0*dwdx37;
    JB[892] = -1.0*dwdx38;
    JB[899] = -1.0*dwdx39;
    JB[907] = 1.0*dwdx38 + 1.0*dwdx39;
    JB[941] = 1.0*dwdx40;
    JB[949] = -1.0*dwdx40;
    JB[980] = 1.0*dwdx41 - 1.0*dwdx42;
    JB[981] = -1.0*dwdx41;
    JB[983] = 1.0*dwdx42;
    JB[997] = 1.0*dwdx45 + 1.0*dwdx46;
    JB[1014] = -1.0*dwdx43;
    JB[1015] = 1.0*dwdx43 + 1.0*dwdx44;
    JB[1016] = -1.0*dwdx44;
    JB[1017] = -1.0*dwdx45;
    JB[1018] = -1.0*dwdx46;
    JB[1058] = -1.0*dwdx49;
    JB[1062] = -1.0*dwdx50;
    JB[1063] = -1.0*dwdx51;
    JB[1065] = -1.0*dwdx47;
    JB[1082] = -1.0*dwdx52;
    JB[1085] = 1.0*dwdx47 - 1.0*dwdx48 + 1.0*dwdx49 + 1.0*dwdx50 + 1.0*dwdx52;
}