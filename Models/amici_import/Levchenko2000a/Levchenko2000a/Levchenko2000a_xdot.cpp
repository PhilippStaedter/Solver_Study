#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Levchenko2000a(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r21 + 1.0*flux_r22 + 1.0*flux_r23 - 1.0*flux_r27 + 1.0*flux_r28 + 1.0*flux_r29;
    xdot[1] = 1.0*flux_r10 + 1.0*flux_r11 - 1.0*flux_r15 + 1.0*flux_r16 + 1.0*flux_r17 - 1.0*flux_r9;
    xdot[2] = -1.0*flux_r0 + 1.0*flux_r1 + 1.0*flux_r2 - 1.0*flux_r252 + 1.0*flux_r253 + 1.0*flux_r254 - 1.0*flux_r255 + 1.0*flux_r256 + 1.0*flux_r257 - 1.0*flux_r258 + 1.0*flux_r259 + 1.0*flux_r260 - 1.0*flux_r261 + 1.0*flux_r262 + 1.0*flux_r263 - 1.0*flux_r264 + 1.0*flux_r265 + 1.0*flux_r266 - 1.0*flux_r267 + 1.0*flux_r268 + 1.0*flux_r269 - 1.0*flux_r270 + 1.0*flux_r271 + 1.0*flux_r272 - 1.0*flux_r273 + 1.0*flux_r274 + 1.0*flux_r275 - 1.0*flux_r276 + 1.0*flux_r277 + 1.0*flux_r278 - 1.0*flux_r279 + 1.0*flux_r280 + 1.0*flux_r281 - 1.0*flux_r282 + 1.0*flux_r283 + 1.0*flux_r284 - 1.0*flux_r285 + 1.0*flux_r286 + 1.0*flux_r287 - 1.0*flux_r288 + 1.0*flux_r289 + 1.0*flux_r290 - 1.0*flux_r291 + 1.0*flux_r292 + 1.0*flux_r293 - 1.0*flux_r294 + 1.0*flux_r295 + 1.0*flux_r296 - 1.0*flux_r297 + 1.0*flux_r298 + 1.0*flux_r299;
    xdot[3] = -1.0*flux_r3 + 1.0*flux_r4 + 1.0*flux_r5;
    xdot[4] = -1.0*flux_r18 + 1.0*flux_r19 + 1.0*flux_r23 - 1.0*flux_r30 + 1.0*flux_r31 - 1.0*flux_r32 + 1.0*flux_r33 - 1.0*flux_r34 + 1.0*flux_r35 - 1.0*flux_r36 + 1.0*flux_r37 - 1.0*flux_r38 + 1.0*flux_r39 - 1.0*flux_r40 + 1.0*flux_r41 - 1.0*flux_r42 + 1.0*flux_r43 - 1.0*flux_r44 + 1.0*flux_r45 - 1.0*flux_r46 + 1.0*flux_r47 - 1.0*flux_r48 + 1.0*flux_r49 - 1.0*flux_r50 + 1.0*flux_r51 - 1.0*flux_r52 + 1.0*flux_r53;
    xdot[5] = 1.0*flux_r20 - 1.0*flux_r21 + 1.0*flux_r22 - 1.0*flux_r24 + 1.0*flux_r25 + 1.0*flux_r29 - 1.0*flux_r54 + 1.0*flux_r55 - 1.0*flux_r56 + 1.0*flux_r57 - 1.0*flux_r58 + 1.0*flux_r59 - 1.0*flux_r60 + 1.0*flux_r61 - 1.0*flux_r62 + 1.0*flux_r63 - 1.0*flux_r64 + 1.0*flux_r65 - 1.0*flux_r66 + 1.0*flux_r67 - 1.0*flux_r68 + 1.0*flux_r69 - 1.0*flux_r70 + 1.0*flux_r71 - 1.0*flux_r72 + 1.0*flux_r73 - 1.0*flux_r74 + 1.0*flux_r75 - 1.0*flux_r76 + 1.0*flux_r77;
    xdot[6] = -1.0*flux_r100 + 1.0*flux_r101 + 1.0*flux_r26 - 1.0*flux_r27 + 1.0*flux_r28 - 1.0*flux_r78 + 1.0*flux_r79 - 1.0*flux_r80 + 1.0*flux_r81 - 1.0*flux_r82 + 1.0*flux_r83 - 1.0*flux_r84 + 1.0*flux_r85 - 1.0*flux_r86 + 1.0*flux_r87 - 1.0*flux_r88 + 1.0*flux_r89 - 1.0*flux_r90 + 1.0*flux_r91 - 1.0*flux_r92 + 1.0*flux_r93 - 1.0*flux_r94 + 1.0*flux_r95 - 1.0*flux_r96 + 1.0*flux_r97 - 1.0*flux_r98 + 1.0*flux_r99;
    xdot[7] = -1.0*flux_r102 + 1.0*flux_r103 - 1.0*flux_r104 + 1.0*flux_r105 - 1.0*flux_r106 + 1.0*flux_r107 + 1.0*flux_r11 - 1.0*flux_r120 + 1.0*flux_r121 - 1.0*flux_r122 + 1.0*flux_r123 - 1.0*flux_r124 + 1.0*flux_r125 - 1.0*flux_r138 + 1.0*flux_r139 - 1.0*flux_r140 + 1.0*flux_r141 - 1.0*flux_r142 + 1.0*flux_r143 - 1.0*flux_r156 + 1.0*flux_r157 - 1.0*flux_r158 + 1.0*flux_r159 - 1.0*flux_r160 + 1.0*flux_r161 - 1.0*flux_r6 + 1.0*flux_r7;
    xdot[8] = 1.0*flux_r10 - 1.0*flux_r108 + 1.0*flux_r109 - 1.0*flux_r110 + 1.0*flux_r111 - 1.0*flux_r112 + 1.0*flux_r113 - 1.0*flux_r12 - 1.0*flux_r126 + 1.0*flux_r127 - 1.0*flux_r128 + 1.0*flux_r129 + 1.0*flux_r13 - 1.0*flux_r130 + 1.0*flux_r131 - 1.0*flux_r144 + 1.0*flux_r145 - 1.0*flux_r146 + 1.0*flux_r147 - 1.0*flux_r148 + 1.0*flux_r149 - 1.0*flux_r162 + 1.0*flux_r163 - 1.0*flux_r164 + 1.0*flux_r165 - 1.0*flux_r166 + 1.0*flux_r167 + 1.0*flux_r17 + 1.0*flux_r8 - 1.0*flux_r9;
    xdot[9] = -1.0*flux_r114 + 1.0*flux_r115 - 1.0*flux_r116 + 1.0*flux_r117 - 1.0*flux_r118 + 1.0*flux_r119 - 1.0*flux_r132 + 1.0*flux_r133 - 1.0*flux_r134 + 1.0*flux_r135 - 1.0*flux_r136 + 1.0*flux_r137 + 1.0*flux_r14 - 1.0*flux_r15 - 1.0*flux_r150 + 1.0*flux_r151 - 1.0*flux_r152 + 1.0*flux_r153 - 1.0*flux_r154 + 1.0*flux_r155 + 1.0*flux_r16 - 1.0*flux_r168 + 1.0*flux_r169 - 1.0*flux_r170 + 1.0*flux_r171 - 1.0*flux_r172 + 1.0*flux_r173 - 1.0*flux_r18 + 1.0*flux_r19 + 1.0*flux_r20 - 1.0*flux_r24 + 1.0*flux_r25 + 1.0*flux_r26;
    xdot[10] = -1.0*flux_r0 + 1.0*flux_r1 - 1.0*flux_r174 + 1.0*flux_r175 - 1.0*flux_r178 + 1.0*flux_r179 - 1.0*flux_r182 + 1.0*flux_r183 - 1.0*flux_r186 + 1.0*flux_r187 - 1.0*flux_r190 + 1.0*flux_r191 - 1.0*flux_r194 + 1.0*flux_r195 - 1.0*flux_r198 + 1.0*flux_r199 - 1.0*flux_r202 + 1.0*flux_r203 - 1.0*flux_r206 + 1.0*flux_r207 - 1.0*flux_r210 + 1.0*flux_r211 - 1.0*flux_r214 + 1.0*flux_r215 - 1.0*flux_r218 + 1.0*flux_r219 - 1.0*flux_r222 + 1.0*flux_r223 - 1.0*flux_r226 + 1.0*flux_r227 - 1.0*flux_r230 + 1.0*flux_r231 - 1.0*flux_r234 + 1.0*flux_r235 + 1.0*flux_r5;
    xdot[11] = -1.0*flux_r12 + 1.0*flux_r13 + 1.0*flux_r14 - 1.0*flux_r176 + 1.0*flux_r177 - 1.0*flux_r180 + 1.0*flux_r181 - 1.0*flux_r184 + 1.0*flux_r185 - 1.0*flux_r188 + 1.0*flux_r189 - 1.0*flux_r192 + 1.0*flux_r193 - 1.0*flux_r196 + 1.0*flux_r197 + 1.0*flux_r2 - 1.0*flux_r200 + 1.0*flux_r201 - 1.0*flux_r204 + 1.0*flux_r205 - 1.0*flux_r208 + 1.0*flux_r209 - 1.0*flux_r212 + 1.0*flux_r213 - 1.0*flux_r216 + 1.0*flux_r217 - 1.0*flux_r220 + 1.0*flux_r221 - 1.0*flux_r224 + 1.0*flux_r225 - 1.0*flux_r228 + 1.0*flux_r229 - 1.0*flux_r232 + 1.0*flux_r233 - 1.0*flux_r236 + 1.0*flux_r237 - 1.0*flux_r3 + 1.0*flux_r4 - 1.0*flux_r6 + 1.0*flux_r7 + 1.0*flux_r8;
    xdot[12] = 1.0*flux_r18 - 1.0*flux_r19 - 1.0*flux_r20;
    xdot[13] = 1.0*flux_r24 - 1.0*flux_r25 - 1.0*flux_r26;
    xdot[14] = 1.0*flux_r6 - 1.0*flux_r7 - 1.0*flux_r8;
    xdot[15] = 1.0*flux_r12 - 1.0*flux_r13 - 1.0*flux_r14;
    xdot[16] = 1.0*flux_r21 - 1.0*flux_r22 - 1.0*flux_r23;
    xdot[17] = 1.0*flux_r27 - 1.0*flux_r28 - 1.0*flux_r29;
    xdot[18] = -1.0*flux_r10 - 1.0*flux_r11 + 1.0*flux_r9;
    xdot[19] = 1.0*flux_r15 - 1.0*flux_r16 - 1.0*flux_r17;
    xdot[20] = 1.0*flux_r0 - 1.0*flux_r1 - 1.0*flux_r2;
    xdot[21] = 1.0*flux_r3 - 1.0*flux_r4 - 1.0*flux_r5;
    xdot[22] = -1.0*flux_r102 + 1.0*flux_r103 - 1.0*flux_r108 + 1.0*flux_r109 - 1.0*flux_r114 + 1.0*flux_r115 - 1.0*flux_r174 + 1.0*flux_r175 - 1.0*flux_r176 + 1.0*flux_r177 - 1.0*flux_r30 + 1.0*flux_r31 - 1.0*flux_r54 + 1.0*flux_r55 - 1.0*flux_r78 + 1.0*flux_r79;
    xdot[23] = -1.0*flux_r104 + 1.0*flux_r105 - 1.0*flux_r110 + 1.0*flux_r111 - 1.0*flux_r116 + 1.0*flux_r117 + 1.0*flux_r174 - 1.0*flux_r175 - 1.0*flux_r252 + 1.0*flux_r253 - 1.0*flux_r32 + 1.0*flux_r33 - 1.0*flux_r56 + 1.0*flux_r57 - 1.0*flux_r80 + 1.0*flux_r81;
    xdot[24] = -1.0*flux_r106 + 1.0*flux_r107 - 1.0*flux_r112 + 1.0*flux_r113 - 1.0*flux_r118 + 1.0*flux_r119 + 1.0*flux_r176 - 1.0*flux_r177 + 1.0*flux_r254 - 1.0*flux_r34 + 1.0*flux_r35 - 1.0*flux_r58 + 1.0*flux_r59 - 1.0*flux_r82 + 1.0*flux_r83;
    xdot[25] = 1.0*flux_r102 - 1.0*flux_r103 - 1.0*flux_r178 + 1.0*flux_r179 - 1.0*flux_r180 + 1.0*flux_r181 - 1.0*flux_r36 + 1.0*flux_r37 - 1.0*flux_r60 + 1.0*flux_r61 - 1.0*flux_r84 + 1.0*flux_r85;
    xdot[26] = 1.0*flux_r104 - 1.0*flux_r105 + 1.0*flux_r178 - 1.0*flux_r179 - 1.0*flux_r255 + 1.0*flux_r256 - 1.0*flux_r38 + 1.0*flux_r39 - 1.0*flux_r62 + 1.0*flux_r63 - 1.0*flux_r86 + 1.0*flux_r87;
    xdot[27] = 1.0*flux_r106 - 1.0*flux_r107 + 1.0*flux_r180 - 1.0*flux_r181 - 1.0*flux_r244 + 1.0*flux_r257 - 1.0*flux_r40 + 1.0*flux_r41 - 1.0*flux_r64 + 1.0*flux_r65 - 1.0*flux_r88 + 1.0*flux_r89;
    xdot[28] = 1.0*flux_r108 - 1.0*flux_r109 - 1.0*flux_r182 + 1.0*flux_r183 - 1.0*flux_r184 + 1.0*flux_r185 - 1.0*flux_r42 + 1.0*flux_r43 - 1.0*flux_r66 + 1.0*flux_r67 - 1.0*flux_r90 + 1.0*flux_r91;
    xdot[29] = 1.0*flux_r110 - 1.0*flux_r111 + 1.0*flux_r182 - 1.0*flux_r183 - 1.0*flux_r258 + 1.0*flux_r259 - 1.0*flux_r44 + 1.0*flux_r45 - 1.0*flux_r68 + 1.0*flux_r69 - 1.0*flux_r92 + 1.0*flux_r93;
    xdot[30] = 1.0*flux_r112 - 1.0*flux_r113 + 1.0*flux_r184 - 1.0*flux_r185 + 1.0*flux_r244 - 1.0*flux_r245 + 1.0*flux_r260 - 1.0*flux_r46 + 1.0*flux_r47 - 1.0*flux_r70 + 1.0*flux_r71 - 1.0*flux_r94 + 1.0*flux_r95;
    xdot[31] = 1.0*flux_r114 - 1.0*flux_r115 - 1.0*flux_r186 + 1.0*flux_r187 - 1.0*flux_r188 + 1.0*flux_r189 - 1.0*flux_r48 + 1.0*flux_r49 - 1.0*flux_r72 + 1.0*flux_r73 - 1.0*flux_r96 + 1.0*flux_r97;
    xdot[32] = 1.0*flux_r116 - 1.0*flux_r117 + 1.0*flux_r186 - 1.0*flux_r187 - 1.0*flux_r261 + 1.0*flux_r262 - 1.0*flux_r50 + 1.0*flux_r51 - 1.0*flux_r74 + 1.0*flux_r75 - 1.0*flux_r98 + 1.0*flux_r99;
    xdot[33] = -1.0*flux_r100 + 1.0*flux_r101 + 1.0*flux_r118 - 1.0*flux_r119 + 1.0*flux_r188 - 1.0*flux_r189 + 1.0*flux_r245 + 1.0*flux_r263 - 1.0*flux_r52 + 1.0*flux_r53 - 1.0*flux_r76 + 1.0*flux_r77;
    xdot[34] = -1.0*flux_r120 + 1.0*flux_r121 - 1.0*flux_r126 + 1.0*flux_r127 - 1.0*flux_r132 + 1.0*flux_r133 - 1.0*flux_r190 + 1.0*flux_r191 - 1.0*flux_r192 + 1.0*flux_r193 + 1.0*flux_r30 - 1.0*flux_r31;
    xdot[35] = -1.0*flux_r122 + 1.0*flux_r123 - 1.0*flux_r128 + 1.0*flux_r129 - 1.0*flux_r134 + 1.0*flux_r135 + 1.0*flux_r190 - 1.0*flux_r191 - 1.0*flux_r264 + 1.0*flux_r265 + 1.0*flux_r32 - 1.0*flux_r33;
    xdot[36] = -1.0*flux_r124 + 1.0*flux_r125 - 1.0*flux_r130 + 1.0*flux_r131 - 1.0*flux_r136 + 1.0*flux_r137 + 1.0*flux_r192 - 1.0*flux_r193 + 1.0*flux_r266 + 1.0*flux_r34 - 1.0*flux_r35;
    xdot[37] = 1.0*flux_r120 - 1.0*flux_r121 - 1.0*flux_r194 + 1.0*flux_r195 - 1.0*flux_r196 + 1.0*flux_r197 + 1.0*flux_r36 - 1.0*flux_r37;
    xdot[38] = 1.0*flux_r122 - 1.0*flux_r123 + 1.0*flux_r194 - 1.0*flux_r195 - 1.0*flux_r267 + 1.0*flux_r268 + 1.0*flux_r38 - 1.0*flux_r39;
    xdot[39] = 1.0*flux_r124 - 1.0*flux_r125 + 1.0*flux_r196 - 1.0*flux_r197 - 1.0*flux_r246 + 1.0*flux_r269 + 1.0*flux_r40 - 1.0*flux_r41;
    xdot[40] = 1.0*flux_r126 - 1.0*flux_r127 - 1.0*flux_r198 + 1.0*flux_r199 - 1.0*flux_r200 + 1.0*flux_r201 + 1.0*flux_r42 - 1.0*flux_r43;
    xdot[41] = 1.0*flux_r128 - 1.0*flux_r129 + 1.0*flux_r198 - 1.0*flux_r199 - 1.0*flux_r270 + 1.0*flux_r271 + 1.0*flux_r44 - 1.0*flux_r45;
    xdot[42] = 1.0*flux_r130 - 1.0*flux_r131 + 1.0*flux_r200 - 1.0*flux_r201 + 1.0*flux_r246 - 1.0*flux_r247 + 1.0*flux_r272 + 1.0*flux_r46 - 1.0*flux_r47;
    xdot[43] = 1.0*flux_r132 - 1.0*flux_r133 - 1.0*flux_r202 + 1.0*flux_r203 - 1.0*flux_r204 + 1.0*flux_r205 - 1.0*flux_r238 + 1.0*flux_r48 - 1.0*flux_r49;
    xdot[44] = 1.0*flux_r134 - 1.0*flux_r135 + 1.0*flux_r202 - 1.0*flux_r203 - 1.0*flux_r239 - 1.0*flux_r273 + 1.0*flux_r274 + 1.0*flux_r50 - 1.0*flux_r51;
    xdot[45] = 1.0*flux_r136 - 1.0*flux_r137 + 1.0*flux_r204 - 1.0*flux_r205 - 1.0*flux_r240 + 1.0*flux_r247 + 1.0*flux_r275 + 1.0*flux_r52 - 1.0*flux_r53;
    xdot[46] = -1.0*flux_r138 + 1.0*flux_r139 - 1.0*flux_r144 + 1.0*flux_r145 - 1.0*flux_r150 + 1.0*flux_r151 - 1.0*flux_r206 + 1.0*flux_r207 - 1.0*flux_r208 + 1.0*flux_r209 + 1.0*flux_r54 - 1.0*flux_r55;
    xdot[47] = -1.0*flux_r140 + 1.0*flux_r141 - 1.0*flux_r146 + 1.0*flux_r147 - 1.0*flux_r152 + 1.0*flux_r153 + 1.0*flux_r206 - 1.0*flux_r207 - 1.0*flux_r276 + 1.0*flux_r277 + 1.0*flux_r56 - 1.0*flux_r57;
    xdot[48] = -1.0*flux_r142 + 1.0*flux_r143 - 1.0*flux_r148 + 1.0*flux_r149 - 1.0*flux_r154 + 1.0*flux_r155 + 1.0*flux_r208 - 1.0*flux_r209 + 1.0*flux_r278 + 1.0*flux_r58 - 1.0*flux_r59;
    xdot[49] = 1.0*flux_r138 - 1.0*flux_r139 - 1.0*flux_r210 + 1.0*flux_r211 - 1.0*flux_r212 + 1.0*flux_r213 + 1.0*flux_r60 - 1.0*flux_r61;
    xdot[50] = 1.0*flux_r140 - 1.0*flux_r141 + 1.0*flux_r210 - 1.0*flux_r211 - 1.0*flux_r279 + 1.0*flux_r280 + 1.0*flux_r62 - 1.0*flux_r63;
    xdot[51] = 1.0*flux_r142 - 1.0*flux_r143 + 1.0*flux_r212 - 1.0*flux_r213 - 1.0*flux_r248 + 1.0*flux_r281 + 1.0*flux_r64 - 1.0*flux_r65;
    xdot[52] = 1.0*flux_r144 - 1.0*flux_r145 - 1.0*flux_r214 + 1.0*flux_r215 - 1.0*flux_r216 + 1.0*flux_r217 + 1.0*flux_r66 - 1.0*flux_r67;
    xdot[53] = 1.0*flux_r146 - 1.0*flux_r147 + 1.0*flux_r214 - 1.0*flux_r215 - 1.0*flux_r282 + 1.0*flux_r283 + 1.0*flux_r68 - 1.0*flux_r69;
    xdot[54] = 1.0*flux_r148 - 1.0*flux_r149 + 1.0*flux_r216 - 1.0*flux_r217 + 1.0*flux_r248 - 1.0*flux_r249 + 1.0*flux_r284 + 1.0*flux_r70 - 1.0*flux_r71;
    xdot[55] = 1.0*flux_r150 - 1.0*flux_r151 - 1.0*flux_r218 + 1.0*flux_r219 - 1.0*flux_r220 + 1.0*flux_r221 + 1.0*flux_r238 - 1.0*flux_r241 + 1.0*flux_r72 - 1.0*flux_r73;
    xdot[56] = 1.0*flux_r152 - 1.0*flux_r153 + 1.0*flux_r218 - 1.0*flux_r219 + 1.0*flux_r239 - 1.0*flux_r242 - 1.0*flux_r285 + 1.0*flux_r286 + 1.0*flux_r74 - 1.0*flux_r75;
    xdot[57] = 1.0*flux_r154 - 1.0*flux_r155 + 1.0*flux_r220 - 1.0*flux_r221 + 1.0*flux_r240 - 1.0*flux_r243 + 1.0*flux_r249 + 1.0*flux_r287 + 1.0*flux_r76 - 1.0*flux_r77;
    xdot[58] = -1.0*flux_r156 + 1.0*flux_r157 - 1.0*flux_r162 + 1.0*flux_r163 - 1.0*flux_r168 + 1.0*flux_r169 - 1.0*flux_r222 + 1.0*flux_r223 - 1.0*flux_r224 + 1.0*flux_r225 + 1.0*flux_r78 - 1.0*flux_r79;
    xdot[59] = -1.0*flux_r158 + 1.0*flux_r159 - 1.0*flux_r164 + 1.0*flux_r165 - 1.0*flux_r170 + 1.0*flux_r171 + 1.0*flux_r222 - 1.0*flux_r223 - 1.0*flux_r288 + 1.0*flux_r289 + 1.0*flux_r80 - 1.0*flux_r81;
    xdot[60] = -1.0*flux_r160 + 1.0*flux_r161 - 1.0*flux_r166 + 1.0*flux_r167 - 1.0*flux_r172 + 1.0*flux_r173 + 1.0*flux_r224 - 1.0*flux_r225 + 1.0*flux_r290 + 1.0*flux_r82 - 1.0*flux_r83;
    xdot[61] = 1.0*flux_r156 - 1.0*flux_r157 - 1.0*flux_r226 + 1.0*flux_r227 - 1.0*flux_r228 + 1.0*flux_r229 + 1.0*flux_r84 - 1.0*flux_r85;
    xdot[62] = 1.0*flux_r158 - 1.0*flux_r159 + 1.0*flux_r226 - 1.0*flux_r227 - 1.0*flux_r291 + 1.0*flux_r292 + 1.0*flux_r86 - 1.0*flux_r87;
    xdot[63] = 1.0*flux_r160 - 1.0*flux_r161 + 1.0*flux_r228 - 1.0*flux_r229 - 1.0*flux_r250 + 1.0*flux_r293 + 1.0*flux_r88 - 1.0*flux_r89;
    xdot[64] = 1.0*flux_r162 - 1.0*flux_r163 - 1.0*flux_r230 + 1.0*flux_r231 - 1.0*flux_r232 + 1.0*flux_r233 + 1.0*flux_r90 - 1.0*flux_r91;
    xdot[65] = 1.0*flux_r164 - 1.0*flux_r165 + 1.0*flux_r230 - 1.0*flux_r231 - 1.0*flux_r294 + 1.0*flux_r295 + 1.0*flux_r92 - 1.0*flux_r93;
    xdot[66] = 1.0*flux_r166 - 1.0*flux_r167 + 1.0*flux_r232 - 1.0*flux_r233 + 1.0*flux_r250 - 1.0*flux_r251 + 1.0*flux_r296 + 1.0*flux_r94 - 1.0*flux_r95;
    xdot[67] = 1.0*flux_r168 - 1.0*flux_r169 - 1.0*flux_r234 + 1.0*flux_r235 - 1.0*flux_r236 + 1.0*flux_r237 + 1.0*flux_r241 + 1.0*flux_r96 - 1.0*flux_r97;
    xdot[68] = 1.0*flux_r170 - 1.0*flux_r171 + 1.0*flux_r234 - 1.0*flux_r235 + 1.0*flux_r242 - 1.0*flux_r297 + 1.0*flux_r298 + 1.0*flux_r98 - 1.0*flux_r99;
    xdot[69] = 1.0*flux_r100 - 1.0*flux_r101 + 1.0*flux_r172 - 1.0*flux_r173 + 1.0*flux_r236 - 1.0*flux_r237 + 1.0*flux_r243 + 1.0*flux_r251 + 1.0*flux_r299;
    xdot[70] = 1.0*flux_r252 - 1.0*flux_r253 - 1.0*flux_r254;
    xdot[71] = 1.0*flux_r255 - 1.0*flux_r256 - 1.0*flux_r257;
    xdot[72] = 1.0*flux_r258 - 1.0*flux_r259 - 1.0*flux_r260;
    xdot[73] = 1.0*flux_r261 - 1.0*flux_r262 - 1.0*flux_r263;
    xdot[74] = 1.0*flux_r264 - 1.0*flux_r265 - 1.0*flux_r266;
    xdot[75] = 1.0*flux_r267 - 1.0*flux_r268 - 1.0*flux_r269;
    xdot[76] = 1.0*flux_r270 - 1.0*flux_r271 - 1.0*flux_r272;
    xdot[77] = 1.0*flux_r273 - 1.0*flux_r274 - 1.0*flux_r275;
    xdot[78] = 1.0*flux_r276 - 1.0*flux_r277 - 1.0*flux_r278;
    xdot[79] = 1.0*flux_r279 - 1.0*flux_r280 - 1.0*flux_r281;
    xdot[80] = 1.0*flux_r282 - 1.0*flux_r283 - 1.0*flux_r284;
    xdot[81] = 1.0*flux_r285 - 1.0*flux_r286 - 1.0*flux_r287;
    xdot[82] = 1.0*flux_r288 - 1.0*flux_r289 - 1.0*flux_r290;
    xdot[83] = 1.0*flux_r291 - 1.0*flux_r292 - 1.0*flux_r293;
    xdot[84] = 1.0*flux_r294 - 1.0*flux_r295 - 1.0*flux_r296;
    xdot[85] = 1.0*flux_r297 - 1.0*flux_r298 - 1.0*flux_r299;
}