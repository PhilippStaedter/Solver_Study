#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_schilling1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[1] = -1.0*dwdx0;
    JB[34] = 1.0*dwdx1;
    JB[35] = -1.0*dwdx1;
    JB[68] = 1.0*dwdx2;
    JB[69] = -1.0*dwdx2;
    JB[102] = 1.0*dwdx3;
    JB[103] = -1.0*dwdx3;
    JB[136] = 1.0*dwdx4;
    JB[137] = -1.0*dwdx4;
    JB[170] = 1.0*dwdx5;
    JB[171] = -1.0*dwdx5;
    JB[204] = 1.0*dwdx6;
    JB[205] = -1.0*dwdx6;
    JB[238] = 1.0*dwdx7;
    JB[249] = -1.0*dwdx7;
    JB[272] = 1.0*dwdx8 + 1.0*dwdx9;
    JB[285] = -1.0*dwdx8 - 1.0*dwdx9;
    JB[306] = 1.0*dwdx10 + 1.0*dwdx11;
    JB[319] = -1.0*dwdx10 - 1.0*dwdx11;
    JB[342] = 1.0*dwdx12;
    JB[354] = -1.0*dwdx12;
    JB[374] = 1.0*dwdx13;
    JB[386] = -1.0*dwdx13;
    JB[408] = 1.0*dwdx14;
    JB[420] = -1.0*dwdx14;
    JB[442] = 1.0*dwdx15;
    JB[454] = -1.0*dwdx15;
    JB[476] = 1.0*dwdx16;
    JB[488] = -1.0*dwdx16;
    JB[510] = 1.0*dwdx17;
    JB[522] = -1.0*dwdx17;
    JB[544] = 1.0*dwdx18;
    JB[547] = -1.0*dwdx18;
    JB[578] = 1.0*dwdx19;
    JB[581] = -1.0*dwdx19;
    JB[605] = -1.0*dwdx21;
    JB[606] = -1.0*dwdx22;
    JB[610] = -1.0*dwdx20;
    JB[612] = 1.0*dwdx20;
    JB[617] = 1.0*dwdx21;
    JB[618] = 1.0*dwdx22;
    JB[627] = -1.0*dwdx23;
    JB[646] = 1.0*dwdx23;
    JB[675] = 1.0*dwdx25;
    JB[677] = -1.0*dwdx24;
    JB[680] = 1.0*dwdx24 + 1.0*dwdx26 + 1.0*dwdx27;
    JB[687] = -1.0*dwdx25;
    JB[688] = -1.0*dwdx26 - 1.0*dwdx27;
    JB[701] = -1.0*dwdx30;
    JB[714] = 1.0*dwdx28 + 1.0*dwdx29 + 1.0*dwdx30;
    JB[722] = -1.0*dwdx28 - 1.0*dwdx29;
    JB[735] = -1.0*dwdx33;
    JB[748] = 1.0*dwdx31 + 1.0*dwdx32 + 1.0*dwdx33;
    JB[756] = -1.0*dwdx31 - 1.0*dwdx32;
    JB[770] = -1.0*dwdx34;
    JB[775] = 1.0*dwdx36;
    JB[776] = 1.0*dwdx35;
    JB[778] = -1.0*dwdx36;
    JB[779] = -1.0*dwdx35;
    JB[782] = 1.0*dwdx34;
    JB[803] = 1.0*dwdx38;
    JB[804] = -1.0*dwdx37;
    JB[815] = -1.0*dwdx38;
    JB[816] = 1.0*dwdx37;
    JB[838] = -1.0*dwdx40;
    JB[850] = 1.0*dwdx39 + 1.0*dwdx40;
    JB[856] = -1.0*dwdx39;
    JB[872] = -1.0*dwdx42;
    JB[884] = 1.0*dwdx41 + 1.0*dwdx42;
    JB[890] = -1.0*dwdx41;
    JB[904] = 1.0*dwdx45;
    JB[905] = 1.0*dwdx44;
    JB[906] = -1.0*dwdx43;
    JB[916] = -1.0*dwdx45 + 1.0*dwdx47;
    JB[917] = -1.0*dwdx44 + 1.0*dwdx46;
    JB[918] = 1.0*dwdx43;
    JB[922] = -1.0*dwdx47;
    JB[923] = -1.0*dwdx46;
    JB[941] = -1.0*dwdx48;
    JB[952] = 1.0*dwdx48;
    JB[977] = 1.0*dwdx50;
    JB[978] = -1.0*dwdx49;
    JB[985] = -1.0*dwdx50;
    JB[986] = 1.0*dwdx49;
    JB[1010] = 1.0*dwdx52;
    JB[1012] = -1.0*dwdx51;
    JB[1018] = -1.0*dwdx52;
    JB[1020] = 1.0*dwdx51;
    JB[1031] = 1.0*dwdx54;
    JB[1032] = 1.0*dwdx55;
    JB[1044] = -1.0*dwdx54 + 1.0*dwdx56;
    JB[1045] = -1.0*dwdx55 + 1.0*dwdx57;
    JB[1048] = -1.0*dwdx53;
    JB[1052] = -1.0*dwdx56;
    JB[1053] = -1.0*dwdx57;
    JB[1054] = 1.0*dwdx53;
    JB[1064] = 1.0*dwdx59;
    JB[1065] = 1.0*dwdx60;
    JB[1077] = -1.0*dwdx59 + 1.0*dwdx61;
    JB[1078] = -1.0*dwdx60 + 1.0*dwdx62;
    JB[1082] = -1.0*dwdx58;
    JB[1085] = -1.0*dwdx61;
    JB[1086] = -1.0*dwdx62;
    JB[1088] = 1.0*dwdx58;
}