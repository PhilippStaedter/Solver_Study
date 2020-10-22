#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_mcauley1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[35] = -1.0*dwdx1;
    J[36] = 1.0*dwdx2 - 1.0*dwdx6 + 1.0*dwdx8;
    J[58] = -1.0*dwdx37;
    J[69] = -1.0*dwdx1;
    J[70] = -1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 + 1.0*dwdx7;
    J[73] = -1.0*dwdx11;
    J[74] = 1.0*dwdx13;
    J[83] = 1.0*dwdx27 + 1.0*dwdx28;
    J[85] = 1.0*dwdx30;
    J[88] = 1.0*dwdx33;
    J[92] = -1.0*dwdx37;
    J[140] = -1.0*dwdx9;
    J[141] = 1.0*dwdx10;
    J[142] = -1.0*dwdx12;
    J[167] = 1.0*dwdx49;
    J[280] = -1.0*dwdx14 - 1.0*dwdx15;
    J[286] = -1.0*dwdx24;
    J[295] = 1.0*dwdx34;
    J[298] = 1.0*dwdx38;
    J[303] = 1.0*dwdx50;
    J[315] = -1.0*dwdx16;
    J[316] = 1.0*dwdx18;
    J[337] = 1.0*dwdx51;
    J[374] = 1.0*dwdx0;
    J[385] = 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21;
    J[403] = -1.0*dwdx45 - 1.0*dwdx46;
    J[405] = 1.0*dwdx47;
    J[417] = 1.0*dwdx16;
    J[450] = 1.0*dwdx15;
    J[455] = -1.0*dwdx22 - 1.0*dwdx23;
    J[456] = 1.0*dwdx24;
    J[458] = -1.0*dwdx29;
    J[519] = -1.0*dwdx17;
    J[523] = 1.0*dwdx23;
    J[525] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    J[526] = 1.0*dwdx29;
    J[527] = -1.0*dwdx30;
    J[533] = 1.0*dwdx35;
    J[536] = 1.0*dwdx39;
    J[580] = 1.0*dwdx3;
    J[595] = -1.0*dwdx31;
    J[596] = 1.0*dwdx32;
    J[663] = 1.0*dwdx31;
    J[682] = 1.0*dwdx4;
    J[685] = 1.0*dwdx11;
    J[686] = -1.0*dwdx13;
    J[700] = -1.0*dwdx33;
    J[716] = 1.0*dwdx5;
    J[783] = 1.0*dwdx1;
    J[784] = 1.0*dwdx6;
    J[805] = -1.0*dwdx34 - 1.0*dwdx35 - 1.0*dwdx36;
    J[806] = 1.0*dwdx37;
    J[808] = -1.0*dwdx38 - 1.0*dwdx39;
    J[809] = -1.0*dwdx40;
    J[980] = -1.0*dwdx41 + 1.0*dwdx42;
    J[981] = 1.0*dwdx43;
    J[983] = 1.0*dwdx52;
    J[1014] = 1.0*dwdx41;
    J[1015] = -1.0*dwdx43 - 1.0*dwdx44;
    J[1049] = 1.0*dwdx44;
    J[1058] = 1.0*dwdx9;
    J[1059] = -1.0*dwdx10;
    J[1060] = 1.0*dwdx12;
    J[1062] = 1.0*dwdx14;
    J[1063] = 1.0*dwdx17;
    J[1065] = 1.0*dwdx20;
    J[1067] = 1.0*dwdx22;
    J[1069] = 1.0*dwdx25 + 1.0*dwdx26;
    J[1077] = 1.0*dwdx36;
    J[1081] = 1.0*dwdx40;
    J[1082] = -1.0*dwdx42;
    J[1083] = 1.0*dwdx45;
    J[1085] = -1.0*dwdx47 + 1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50 - 1.0*dwdx52;
    J[1099] = 1.0*dwdx21;
    J[1117] = 1.0*dwdx46;
}