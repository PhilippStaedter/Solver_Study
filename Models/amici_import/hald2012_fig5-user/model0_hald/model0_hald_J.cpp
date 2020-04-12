#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_hald(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx1 - 59.0*dwdx2;
    J[1] = -59.0*dwdx5;
    J[8] = 1.0*dwdx24;
    J[12] = 1.0*dwdx31;
    J[20] = 1.0*dwdx48 - 1.0*dwdx49;
    J[21] = 1.0*dwdx52;
    J[25] = 1.0*dwdx59;
    J[31] = 1.0*dwdx2;
    J[32] = -1.0*dwdx3 - 1.0*dwdx4 + 1.0*dwdx5;
    J[50] = -1.0*dwdx46;
    J[61] = -1.0*dwdx63;
    J[64] = 2.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8;
    J[65] = 1.0*dwdx10 + 2.0*dwdx9;
    J[66] = 2.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx15 - 1.0*dwdx16;
    J[69] = -1.0*dwdx23;
    J[72] = 1.0*dwdx27;
    J[74] = 1.0*dwdx33;
    J[76] = 1.0*dwdx38;
    J[86] = -1.0*dwdx57 - 1.0*dwdx58;
    J[95] = -1.0*dwdx6;
    J[96] = -1.0*dwdx9;
    J[97] = -1.0*dwdx11;
    J[126] = -1.0*dwdx6 + 1.0*dwdx7 + 1.0*dwdx8;
    J[127] = -1.0*dwdx10 - 1.0*dwdx9;
    J[128] = -1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15 + 1.0*dwdx16;
    J[131] = 1.0*dwdx23;
    J[134] = -1.0*dwdx27;
    J[136] = -1.0*dwdx33;
    J[138] = -1.0*dwdx38;
    J[148] = 1.0*dwdx57 + 1.0*dwdx58;
    J[160] = 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20;
    J[161] = -1.0*dwdx21;
    J[166] = 1.0*dwdx29;
    J[168] = 1.0*dwdx34 - 1.0*dwdx36;
    J[173] = -1.0*dwdx42;
    J[175] = -1.0*dwdx51;
    J[176] = -1.0*dwdx54;
    J[191] = 1.0*dwdx18;
    J[192] = 1.0*dwdx21;
    J[204] = 1.0*dwdx42;
    J[219] = -1.0*dwdx8;
    J[221] = -1.0*dwdx16;
    J[224] = 1.0*dwdx22 - 1.0*dwdx23;
    J[230] = 1.0*dwdx35;
    J[237] = 1.0*dwdx50;
    J[238] = 1.0*dwdx53;
    J[241] = -1.0*dwdx58;
    J[248] = -1.0*dwdx0;
    J[256] = -1.0*dwdx24 - 59.0*dwdx25;
    J[257] = -59.0*dwdx26;
    J[268] = -1.0*dwdx48;
    J[269] = -1.0*dwdx52;
    J[287] = 1.0*dwdx25;
    J[288] = 1.0*dwdx26;
    J[313] = -1.0*dwdx10;
    J[314] = -1.0*dwdx14;
    J[320] = -1.0*dwdx27 + 1.0*dwdx28;
    J[322] = 1.0*dwdx32;
    J[344] = 1.0*dwdx10;
    J[345] = 1.0*dwdx14;
    J[346] = -1.0*dwdx17;
    J[351] = 1.0*dwdx27;
    J[352] = -1.0*dwdx29;
    J[354] = -1.0*dwdx34;
    J[376] = 1.0*dwdx13 - 1.0*dwdx15;
    J[382] = -1.0*dwdx28;
    J[384] = -1.0*dwdx32 - 1.0*dwdx33;
    J[386] = 1.0*dwdx38;
    J[408] = 1.0*dwdx17 + 1.0*dwdx19;
    J[410] = -1.0*dwdx22;
    J[414] = 1.0*dwdx29;
    J[416] = 1.0*dwdx34 - 1.0*dwdx35 + 1.0*dwdx36;
    J[423] = -1.0*dwdx50;
    J[424] = -1.0*dwdx53;
    J[438] = -1.0*dwdx13;
    J[446] = 59.0*dwdx30;
    J[448] = 59.0*dwdx37 - 1.0*dwdx38;
    J[449] = 59.0*dwdx39;
    J[477] = -1.0*dwdx30;
    J[479] = -1.0*dwdx37;
    J[480] = -1.0*dwdx39;
    J[501] = 1.0*dwdx20;
    J[512] = -59.0*dwdx40;
    J[513] = -59.0*dwdx41;
    J[516] = 1.0*dwdx51;
    J[517] = 1.0*dwdx54;
    J[543] = 1.0*dwdx40;
    J[544] = 1.0*dwdx41;
    J[563] = -1.0*dwdx18;
    J[564] = -1.0*dwdx21;
    J[576] = -1.0*dwdx42 - 1.0*dwdx43 - 59.0*dwdx44;
    J[577] = -59.0*dwdx47;
    J[583] = -1.0*dwdx60;
    J[584] = -1.0*dwdx62;
    J[590] = -1.0*dwdx4;
    J[607] = 1.0*dwdx44;
    J[608] = -1.0*dwdx45 - 1.0*dwdx46 + 1.0*dwdx47;
    J[619] = -1.0*dwdx63;
    J[620] = -1.0*dwdx0 - 1.0*dwdx1;
    J[625] = 1.0*dwdx20;
    J[627] = -1.0*dwdx22;
    J[628] = -1.0*dwdx24;
    J[633] = -1.0*dwdx35;
    J[640] = -1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50 + 1.0*dwdx51;
    J[641] = -1.0*dwdx52 - 1.0*dwdx53 + 1.0*dwdx54;
    J[651] = 1.0*dwdx0 + 1.0*dwdx1;
    J[656] = -1.0*dwdx20;
    J[658] = 1.0*dwdx22;
    J[659] = 1.0*dwdx24;
    J[664] = 1.0*dwdx35;
    J[671] = 1.0*dwdx48 + 1.0*dwdx49 + 1.0*dwdx50 - 1.0*dwdx51;
    J[672] = 1.0*dwdx52 + 1.0*dwdx53 - 1.0*dwdx54;
    J[682] = 1.0*dwdx1;
    J[702] = 1.0*dwdx49;
    J[704] = -59.0*dwdx55;
    J[705] = -59.0*dwdx56;
    J[735] = 1.0*dwdx55;
    J[736] = 1.0*dwdx56;
    J[746] = -1.0*dwdx7 + 1.0*dwdx8;
    J[748] = 1.0*dwdx16;
    J[751] = 1.0*dwdx23;
    J[768] = -1.0*dwdx57 + 1.0*dwdx58;
    J[777] = 1.0*dwdx7;
    J[787] = -1.0*dwdx31;
    J[793] = -1.0*dwdx43;
    J[799] = 1.0*dwdx57;
    J[800] = -1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61;
    J[801] = -1.0*dwdx62;
    J[824] = 1.0*dwdx43;
    J[831] = 1.0*dwdx60;
    J[832] = 1.0*dwdx62;
    J[903] = 1.0*dwdx15;
    J[911] = 1.0*dwdx33;
    J[931] = 1.0*dwdx4;
    J[949] = 1.0*dwdx46;
    J[960] = 1.0*dwdx63;
}