#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_Ouzounoglou2014(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8;
    J[3] = 1.0*dwdx11;
    J[4] = 1.0*dwdx13;
    J[15] = 1.0*dwdx22;
    J[16] = 1.0*dwdx24;
    J[17] = 1.0*dwdx26;
    J[18] = 1.0*dwdx28;
    J[19] = 1.0*dwdx30;
    J[20] = 1.0*dwdx32 + 1.0*dwdx33;
    J[21] = 1.0*dwdx35;
    J[57] = -1.0*dwdx86;
    J[58] = -1.0*dwdx90;
    J[59] = -1.0*dwdx111;
    J[60] = -1.0*dwdx116;
    J[63] = -1.0*dwdx127;
    J[64] = -1.0*dwdx132;
    J[68] = -1.0*dwdx145;
    J[69] = 1.0*dwdx151;
    J[70] = -1.0*dwdx156;
    J[71] = -1.0*dwdx161;
    J[91] = -1.0*dwdx9;
    J[93] = 1.0*dwdx11;
    J[182] = -1.0*dwdx10;
    J[184] = 1.0*dwdx13;
    J[270] = 1.0*dwdx0;
    J[273] = -1.0*dwdx11 - 1.0*dwdx12;
    J[328] = 1.0*dwdx90 - 1.0*dwdx99;
    J[360] = 1.0*dwdx1;
    J[363] = 1.0*dwdx12;
    J[364] = -1.0*dwdx13 - 1.0*dwdx14;
    J[418] = -1.0*dwdx100 + 1.0*dwdx99;
    J[419] = 1.0*dwdx111;
    J[452] = 1.0*dwdx10;
    J[541] = 1.0*dwdx9;
    J[637] = -1.0*dwdx15;
    J[652] = 1.0*dwdx36;
    J[687] = 1.0*dwdx71 - 1.0*dwdx82;
    J[727] = 1.0*dwdx15;
    J[728] = -1.0*dwdx16;
    J[777] = 1.0*dwdx82 - 1.0*dwdx83;
    J[818] = 1.0*dwdx16;
    J[819] = -1.0*dwdx17;
    J[867] = 1.0*dwdx83 - 1.0*dwdx84;
    J[910] = -1.0*dwdx18;
    J[911] = 1.0*dwdx19;
    J[957] = 1.0*dwdx72 - 1.0*dwdx73;
    J[1001] = -1.0*dwdx19;
    J[1002] = 1.0*dwdx20;
    J[1047] = -1.0*dwdx72 + 1.0*dwdx85;
    J[1089] = 1.0*dwdx17;
    J[1092] = -1.0*dwdx20;
    J[1137] = 1.0*dwdx84 - 1.0*dwdx85;
    J[1180] = 1.0*dwdx18;
    J[1183] = -1.0*dwdx21;
    J[1227] = 1.0*dwdx73 - 1.0*dwdx74;
    J[1273] = 1.0*dwdx21;
    J[1317] = 1.0*dwdx74;
    J[1350] = 1.0*dwdx5;
    J[1354] = 1.0*dwdx14;
    J[1365] = -1.0*dwdx22 - 1.0*dwdx23;
    J[1408] = 1.0*dwdx100 - 1.0*dwdx101;
    J[1410] = 1.0*dwdx116;
    J[1440] = 1.0*dwdx6;
    J[1455] = 1.0*dwdx23;
    J[1456] = -1.0*dwdx24 - 1.0*dwdx25;
    J[1498] = 1.0*dwdx101 - 1.0*dwdx102;
    J[1504] = 1.0*dwdx132;
    J[1530] = 1.0*dwdx2;
    J[1546] = 1.0*dwdx25;
    J[1547] = -1.0*dwdx26 - 1.0*dwdx27;
    J[1588] = 1.0*dwdx102 - 1.0*dwdx103;
    J[1593] = 1.0*dwdx127;
    J[1620] = 1.0*dwdx7;
    J[1637] = 1.0*dwdx27;
    J[1638] = -1.0*dwdx28 - 1.0*dwdx29;
    J[1678] = 1.0*dwdx103 - 1.0*dwdx104;
    J[1691] = 1.0*dwdx161;
    J[1710] = 1.0*dwdx3;
    J[1728] = 1.0*dwdx29;
    J[1729] = -1.0*dwdx30 - 1.0*dwdx31;
    J[1768] = 1.0*dwdx104 - 1.0*dwdx105;
    J[1780] = 1.0*dwdx156;
    J[1819] = 1.0*dwdx31;
    J[1820] = -1.0*dwdx32 - 1.0*dwdx33 - 1.0*dwdx34;
    J[1858] = 1.0*dwdx105 - 1.0*dwdx106;
    J[1869] = -1.0*dwdx151;
    J[1890] = 1.0*dwdx4;
    J[1910] = 1.0*dwdx34;
    J[1911] = -1.0*dwdx35;
    J[1948] = 1.0*dwdx106;
    J[1958] = 1.0*dwdx145;
    J[1980] = 1.0*dwdx8;
    J[2002] = -1.0*dwdx36;
    J[2037] = -1.0*dwdx71 + 1.0*dwdx86;
    J[2109] = 1.0*dwdx41;
    J[2119] = 1.0*dwdx50;
    J[2198] = 1.0*dwdx40;
    J[2285] = 1.0*dwdx37;
    J[2376] = 1.0*dwdx38;
    J[2467] = 1.0*dwdx39;
    J[2560] = 1.0*dwdx42;
    J[2651] = 1.0*dwdx43;
    J[2750] = 1.0*dwdx51;
    J[2751] = 1.0*dwdx52;
    J[2838] = 1.0*dwdx49;
    J[2927] = 1.0*dwdx48;
    J[3016] = 1.0*dwdx47;
    J[3103] = 1.0*dwdx44;
    J[3185] = -1.0*dwdx37;
    J[3214] = 1.0*dwdx129;
    J[3276] = -1.0*dwdx38;
    J[3303] = 1.0*dwdx124;
    J[3367] = -1.0*dwdx39;
    J[3401] = 1.0*dwdx158;
    J[3458] = -1.0*dwdx40;
    J[3480] = 1.0*dwdx113;
    J[3549] = -1.0*dwdx41;
    J[3569] = 1.0*dwdx109;
    J[3640] = -1.0*dwdx42;
    J[3670] = 1.0*dwdx153;
    J[3731] = -1.0*dwdx43;
    J[3759] = 1.0*dwdx148;
    J[3824] = 1.0*dwdx45;
    J[3825] = 1.0*dwdx46;
    J[3913] = -1.0*dwdx44;
    J[3936] = 1.0*dwdx138;
    J[4004] = -1.0*dwdx45;
    J[4025] = 1.0*dwdx134;
    J[4095] = -1.0*dwdx46;
    J[4111] = 1.0*dwdx118;
    J[4186] = -1.0*dwdx47;
    J[4192] = 1.0*dwdx53;
    J[4277] = -1.0*dwdx48;
    J[4283] = 1.0*dwdx57;
    J[4368] = -1.0*dwdx49;
    J[4375] = 1.0*dwdx63;
    J[4459] = -1.0*dwdx50;
    J[4468] = 1.0*dwdx108;
    J[4550] = -1.0*dwdx51;
    J[4557] = 1.0*dwdx87;
    J[4641] = -1.0*dwdx52;
    J[4646] = 1.0*dwdx67;
    J[4732] = -1.0*dwdx53 - 1.0*dwdx54 - 1.0*dwdx55 - 1.0*dwdx56;
    J[4733] = 1.0*dwdx58;
    J[4737] = 1.0*dwdx77 - 1.0*dwdx78;
    J[4741] = 1.0*dwdx121;
    J[4753] = -1.0*dwdx174;
    J[4822] = 1.0*dwdx56;
    J[4823] = -1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60;
    J[4825] = 1.0*dwdx64;
    J[4827] = 1.0*dwdx76 - 1.0*dwdx77;
    J[4843] = -1.0*dwdx173;
    J[5003] = 1.0*dwdx60;
    J[5005] = -1.0*dwdx63 - 1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66;
    J[5006] = 1.0*dwdx68;
    J[5007] = 1.0*dwdx75 - 1.0*dwdx76;
    J[5023] = -1.0*dwdx172;
    J[5095] = 1.0*dwdx66;
    J[5096] = -1.0*dwdx67 - 1.0*dwdx68 - 1.0*dwdx69;
    J[5097] = 1.0*dwdx70 - 1.0*dwdx75;
    J[5130] = -1.0*dwdx8;
    J[5137] = -1.0*dwdx15;
    J[5138] = -1.0*dwdx16;
    J[5139] = -1.0*dwdx17;
    J[5140] = -1.0*dwdx18;
    J[5141] = -1.0*dwdx19;
    J[5142] = -1.0*dwdx20;
    J[5143] = -1.0*dwdx21;
    J[5152] = -1.0*dwdx36;
    J[5182] = -1.0*dwdx54 + 1.0*dwdx56;
    J[5183] = -1.0*dwdx58 + 1.0*dwdx60;
    J[5185] = -1.0*dwdx64 + 1.0*dwdx66;
    J[5186] = -1.0*dwdx68 + 2.0*dwdx69;
    J[5187] = -2.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79 - 1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx82 - 1.0*dwdx83 - 1.0*dwdx84 - 1.0*dwdx85 - 1.0*dwdx86 - 1.0*dwdx87;
    J[5188] = 1.0*dwdx89;
    J[5191] = -1.0*dwdx119 + 1.0*dwdx121;
    J[5192] = 1.0*dwdx122;
    J[5195] = -1.0*dwdx135 + 1.0*dwdx137;
    J[5196] = -1.0*dwdx139 + 1.0*dwdx141;
    J[5197] = 1.0*dwdx143;
    J[5220] = -1.0*dwdx0;
    J[5223] = -1.0*dwdx12;
    J[5224] = -1.0*dwdx14;
    J[5235] = -1.0*dwdx23;
    J[5236] = -1.0*dwdx25;
    J[5237] = -1.0*dwdx27;
    J[5238] = -1.0*dwdx29;
    J[5239] = -1.0*dwdx31;
    J[5240] = -1.0*dwdx34;
    J[5274] = 1.0*dwdx61;
    J[5278] = -1.0*dwdx100 - 1.0*dwdx101 - 1.0*dwdx102 - 1.0*dwdx103 - 1.0*dwdx104 - 1.0*dwdx105 - 1.0*dwdx106 - 1.0*dwdx107 - 1.0*dwdx108 - 2.0*dwdx88 - 1.0*dwdx89 - 1.0*dwdx90 - 1.0*dwdx91 - 1.0*dwdx92 - 1.0*dwdx93 - 1.0*dwdx94 - 1.0*dwdx95 - 1.0*dwdx96 - 1.0*dwdx97 - 1.0*dwdx98 - 1.0*dwdx99;
    J[5279] = -1.0*dwdx110 + 2.0*dwdx112;
    J[5280] = -1.0*dwdx114 + 1.0*dwdx117;
    J[5282] = -1.0*dwdx122;
    J[5283] = -1.0*dwdx125 + 1.0*dwdx128;
    J[5284] = -1.0*dwdx130 + 1.0*dwdx133;
    J[5288] = -1.0*dwdx146 + 1.0*dwdx147;
    J[5289] = -1.0*dwdx149 + 1.0*dwdx152;
    J[5290] = -1.0*dwdx154 + 1.0*dwdx157;
    J[5291] = -1.0*dwdx159 + 1.0*dwdx162;
    J[5292] = -1.0*dwdx164;
    J[5310] = -1.0*dwdx1;
    J[5368] = 1.0*dwdx88 - 1.0*dwdx91;
    J[5369] = -1.0*dwdx109 - 1.0*dwdx110 - 1.0*dwdx111 - 1.0*dwdx112;
    J[5370] = 1.0*dwdx117;
    J[5400] = -1.0*dwdx5;
    J[5415] = 1.0*dwdx22;
    J[5458] = 1.0*dwdx91 - 1.0*dwdx92;
    J[5459] = 1.0*dwdx110;
    J[5460] = -1.0*dwdx113 - 1.0*dwdx114 - 1.0*dwdx115 - 1.0*dwdx116 - 1.0*dwdx117;
    J[5464] = 1.0*dwdx133;
    J[5473] = -1.0*dwdx165;
    J[5542] = 1.0*dwdx54;
    J[5547] = 1.0*dwdx78 - 1.0*dwdx79;
    J[5551] = -1.0*dwdx118 - 1.0*dwdx119 - 1.0*dwdx120 - 1.0*dwdx121;
    J[5555] = 1.0*dwdx137;
    J[5563] = -1.0*dwdx175;
    J[5634] = 1.0*dwdx62;
    J[5638] = -1.0*dwdx89;
    J[5642] = -1.0*dwdx122 - 1.0*dwdx123;
    J[5670] = -1.0*dwdx2;
    J[5687] = 1.0*dwdx26;
    J[5728] = 1.0*dwdx93 - 1.0*dwdx94;
    J[5733] = -1.0*dwdx124 - 1.0*dwdx125 - 1.0*dwdx126 - 1.0*dwdx127 - 1.0*dwdx128;
    J[5734] = 1.0*dwdx130;
    J[5741] = 1.0*dwdx162;
    J[5743] = -1.0*dwdx167;
    J[5760] = -1.0*dwdx6;
    J[5776] = 1.0*dwdx24;
    J[5818] = 1.0*dwdx92 - 1.0*dwdx93;
    J[5820] = 1.0*dwdx114;
    J[5823] = 1.0*dwdx128;
    J[5824] = -1.0*dwdx129 - 1.0*dwdx130 - 1.0*dwdx131 - 1.0*dwdx132 - 1.0*dwdx133;
    J[5833] = -1.0*dwdx166;
    J[5907] = 1.0*dwdx79 - 1.0*dwdx80;
    J[5911] = 1.0*dwdx119;
    J[5915] = -1.0*dwdx134 - 1.0*dwdx135 - 1.0*dwdx136 - 1.0*dwdx137;
    J[5916] = 1.0*dwdx141;
    J[5923] = -1.0*dwdx176;
    J[5997] = 1.0*dwdx80 - 1.0*dwdx81;
    J[6005] = 1.0*dwdx135;
    J[6006] = -1.0*dwdx138 - 1.0*dwdx139 - 1.0*dwdx140 - 1.0*dwdx141;
    J[6007] = 1.0*dwdx143;
    J[6013] = -1.0*dwdx177;
    J[6087] = 1.0*dwdx81;
    J[6096] = 1.0*dwdx139;
    J[6097] = -1.0*dwdx142 - 1.0*dwdx143;
    J[6103] = -1.0*dwdx178;
    J[6120] = -1.0*dwdx4;
    J[6141] = 1.0*dwdx35;
    J[6178] = 1.0*dwdx97 - 1.0*dwdx98;
    J[6188] = -1.0*dwdx144 - 1.0*dwdx145 - 1.0*dwdx146 - 1.0*dwdx147;
    J[6189] = 1.0*dwdx149;
    J[6193] = -1.0*dwdx171;
    J[6230] = 1.0*dwdx32 - 1.0*dwdx33;
    J[6268] = 1.0*dwdx96 - 1.0*dwdx97;
    J[6278] = 1.0*dwdx147;
    J[6279] = -1.0*dwdx148 - 1.0*dwdx149 - 1.0*dwdx150 - 1.0*dwdx151 - 1.0*dwdx152;
    J[6280] = 1.0*dwdx154;
    J[6283] = -1.0*dwdx170;
    J[6300] = -1.0*dwdx3;
    J[6319] = 1.0*dwdx30;
    J[6358] = 1.0*dwdx95 - 1.0*dwdx96;
    J[6369] = 1.0*dwdx152;
    J[6370] = -1.0*dwdx153 - 1.0*dwdx154 - 1.0*dwdx155 - 1.0*dwdx156 - 1.0*dwdx157;
    J[6371] = 1.0*dwdx159;
    J[6373] = -1.0*dwdx169;
    J[6390] = -1.0*dwdx7;
    J[6408] = 1.0*dwdx28;
    J[6448] = 1.0*dwdx94 - 1.0*dwdx95;
    J[6453] = 1.0*dwdx125;
    J[6460] = 1.0*dwdx157;
    J[6461] = -1.0*dwdx158 - 1.0*dwdx159 - 1.0*dwdx160 - 1.0*dwdx161 - 1.0*dwdx162;
    J[6463] = -1.0*dwdx168;
    J[6538] = 1.0*dwdx98;
    J[6548] = 1.0*dwdx146;
    J[6552] = -1.0*dwdx163;
    J[6553] = -1.0*dwdx179;
    J[6622] = -1.0*dwdx55;
    J[6623] = -1.0*dwdx59;
    J[6625] = -1.0*dwdx65;
    J[6630] = -1.0*dwdx115;
    J[6631] = -1.0*dwdx120;
    J[6633] = -1.0*dwdx126;
    J[6634] = -1.0*dwdx131;
    J[6635] = -1.0*dwdx136;
    J[6636] = -1.0*dwdx140;
    J[6637] = -1.0*dwdx142;
    J[6638] = -1.0*dwdx144;
    J[6639] = -1.0*dwdx150;
    J[6640] = -1.0*dwdx155;
    J[6641] = -1.0*dwdx160;
    J[6642] = -1.0*dwdx163;
    J[6643] = -1.0*dwdx165 - 1.0*dwdx166 - 1.0*dwdx167 - 1.0*dwdx168 - 1.0*dwdx169 - 1.0*dwdx170 - 1.0*dwdx171 - 1.0*dwdx172 - 1.0*dwdx173 - 1.0*dwdx174 - 1.0*dwdx175 - 1.0*dwdx176 - 1.0*dwdx177 - 1.0*dwdx178 - 1.0*dwdx179;
    J[6645] = 1.0*dwdx180;
    J[6646] = 1.0*dwdx181;
    J[6647] = 1.0*dwdx182;
    J[6648] = 1.0*dwdx183;
    J[6649] = 1.0*dwdx184;
    J[6650] = 1.0*dwdx185;
    J[6651] = 1.0*dwdx186;
    J[6652] = 1.0*dwdx187;
    J[6653] = 1.0*dwdx188;
    J[6654] = 1.0*dwdx189;
    J[6655] = 1.0*dwdx190;
    J[6656] = 1.0*dwdx191;
    J[6657] = 1.0*dwdx192;
    J[6658] = 1.0*dwdx193;
    J[6722] = 1.0*dwdx123;
    J[6810] = 1.0*dwdx115;
    J[6823] = 1.0*dwdx165;
    J[6825] = -1.0*dwdx180;
    J[6904] = 1.0*dwdx131;
    J[6913] = 1.0*dwdx166;
    J[6916] = -1.0*dwdx181;
    J[6993] = 1.0*dwdx126;
    J[7003] = 1.0*dwdx167;
    J[7007] = -1.0*dwdx182;
    J[7091] = 1.0*dwdx160;
    J[7093] = 1.0*dwdx168;
    J[7098] = -1.0*dwdx183;
    J[7180] = 1.0*dwdx155;
    J[7183] = 1.0*dwdx169;
    J[7189] = -1.0*dwdx184;
    J[7269] = 1.0*dwdx150;
    J[7273] = 1.0*dwdx170;
    J[7280] = -1.0*dwdx185;
    J[7358] = 1.0*dwdx144;
    J[7363] = 1.0*dwdx171;
    J[7371] = -1.0*dwdx186;
    J[7435] = 1.0*dwdx65;
    J[7453] = 1.0*dwdx172;
    J[7462] = -1.0*dwdx187;
    J[7523] = 1.0*dwdx59;
    J[7543] = 1.0*dwdx173;
    J[7553] = -1.0*dwdx188;
    J[7612] = 1.0*dwdx55;
    J[7633] = 1.0*dwdx174;
    J[7644] = -1.0*dwdx189;
    J[7711] = 1.0*dwdx120;
    J[7723] = 1.0*dwdx175;
    J[7735] = -1.0*dwdx190;
    J[7805] = 1.0*dwdx136;
    J[7813] = 1.0*dwdx176;
    J[7826] = -1.0*dwdx191;
    J[7896] = 1.0*dwdx140;
    J[7903] = 1.0*dwdx177;
    J[7917] = -1.0*dwdx192;
    J[7987] = 1.0*dwdx142;
    J[7993] = 1.0*dwdx178;
    J[8008] = -1.0*dwdx193;
    J[8082] = 1.0*dwdx163;
    J[8083] = 1.0*dwdx179;
}