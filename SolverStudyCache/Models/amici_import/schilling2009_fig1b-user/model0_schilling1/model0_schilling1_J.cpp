#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_schilling1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[19] = 1.0*dwdx23;
    J[33] = 1.0*dwdx0;
    J[34] = -1.0*dwdx1;
    J[67] = 1.0*dwdx1;
    J[68] = -1.0*dwdx2;
    J[101] = 1.0*dwdx2;
    J[102] = -1.0*dwdx3;
    J[135] = 1.0*dwdx3;
    J[136] = -1.0*dwdx4;
    J[169] = 1.0*dwdx4;
    J[170] = -1.0*dwdx5;
    J[203] = 1.0*dwdx5;
    J[204] = -1.0*dwdx6;
    J[237] = 1.0*dwdx6;
    J[238] = -1.0*dwdx7;
    J[272] = -1.0*dwdx8 - 1.0*dwdx9;
    J[285] = 1.0*dwdx30;
    J[295] = -1.0*dwdx54;
    J[296] = -1.0*dwdx59;
    J[306] = -1.0*dwdx10 - 1.0*dwdx11;
    J[319] = 1.0*dwdx33;
    J[328] = -1.0*dwdx55;
    J[329] = -1.0*dwdx60;
    J[374] = -1.0*dwdx13;
    J[381] = 1.0*dwdx21;
    J[386] = 1.0*dwdx34;
    J[387] = -1.0*dwdx38;
    J[406] = -1.0*dwdx12;
    J[408] = -1.0*dwdx14;
    J[414] = 1.0*dwdx22;
    J[420] = 1.0*dwdx37;
    J[442] = -1.0*dwdx15;
    J[454] = 1.0*dwdx40;
    J[456] = -1.0*dwdx45;
    J[476] = -1.0*dwdx16;
    J[488] = 1.0*dwdx42;
    J[489] = -1.0*dwdx44;
    J[510] = -1.0*dwdx17;
    J[515] = -1.0*dwdx25;
    J[522] = 1.0*dwdx43;
    J[544] = -1.0*dwdx18;
    J[546] = 1.0*dwdx20;
    J[551] = -1.0*dwdx36;
    J[578] = -1.0*dwdx19;
    J[581] = 1.0*dwdx24;
    J[584] = -1.0*dwdx35;
    J[589] = 1.0*dwdx48;
    J[601] = 1.0*dwdx7;
    J[612] = -1.0*dwdx20;
    J[643] = 1.0*dwdx18;
    J[646] = -1.0*dwdx23;
    J[650] = 1.0*dwdx36;
    J[677] = 1.0*dwdx19;
    J[680] = -1.0*dwdx24 - 1.0*dwdx26 - 1.0*dwdx27;
    J[683] = 1.0*dwdx35;
    J[689] = -1.0*dwdx50;
    J[690] = -1.0*dwdx52;
    J[701] = 1.0*dwdx8 + 1.0*dwdx9;
    J[714] = -1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30;
    J[722] = 1.0*dwdx49;
    J[724] = 1.0*dwdx54 - 1.0*dwdx56;
    J[725] = 1.0*dwdx59 - 1.0*dwdx61;
    J[735] = 1.0*dwdx10 + 1.0*dwdx11;
    J[748] = -1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33;
    J[756] = 1.0*dwdx51;
    J[757] = 1.0*dwdx55 - 1.0*dwdx57;
    J[758] = 1.0*dwdx60 - 1.0*dwdx62;
    J[770] = 1.0*dwdx13;
    J[777] = -1.0*dwdx21;
    J[782] = -1.0*dwdx34;
    J[783] = 1.0*dwdx38;
    J[802] = 1.0*dwdx12;
    J[804] = 1.0*dwdx14;
    J[810] = -1.0*dwdx22;
    J[816] = -1.0*dwdx37;
    J[838] = 1.0*dwdx15;
    J[850] = -1.0*dwdx39 - 1.0*dwdx40;
    J[852] = 1.0*dwdx45 - 1.0*dwdx47;
    J[856] = 1.0*dwdx53;
    J[872] = 1.0*dwdx16;
    J[884] = -1.0*dwdx41 - 1.0*dwdx42;
    J[885] = 1.0*dwdx44 - 1.0*dwdx46;
    J[890] = 1.0*dwdx58;
    J[906] = 1.0*dwdx17;
    J[911] = 1.0*dwdx25;
    J[918] = -1.0*dwdx43;
    J[944] = 1.0*dwdx26 + 1.0*dwdx27;
    J[952] = -1.0*dwdx48;
    J[953] = 1.0*dwdx50;
    J[954] = 1.0*dwdx52;
    J[978] = 1.0*dwdx28 + 1.0*dwdx29;
    J[986] = -1.0*dwdx49;
    J[988] = 1.0*dwdx56;
    J[989] = 1.0*dwdx61;
    J[1012] = 1.0*dwdx31 + 1.0*dwdx32;
    J[1020] = -1.0*dwdx51;
    J[1021] = 1.0*dwdx57;
    J[1022] = 1.0*dwdx62;
    J[1048] = 1.0*dwdx39;
    J[1050] = 1.0*dwdx47;
    J[1054] = -1.0*dwdx53;
    J[1082] = 1.0*dwdx41;
    J[1083] = 1.0*dwdx46;
    J[1088] = -1.0*dwdx58;
}