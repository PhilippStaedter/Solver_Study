#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Pathak2013(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*kass_re1*s1 - 1.0*kdiss_re1*s7;
    w[1] = 1.0*kass_re2*s2 - 1.0*kdiss_re2*s7;
    w[2] = 1.0*kass_re3*s2 - 1.0*kdiss_re3*s8;
    w[3] = 1.0*kass_re4*s2 - 1.0*kdiss_re4*s9;
    w[4] = 1.0*kass_re5*s6 - 1.0*kdiss_re5*s12;
    w[5] = 1.0*kass_re6*s6 - 1.0*kdiss_re6*s11;
    w[6] = 1.0*kass_re7*s6 - 1.0*kdiss_re7*s10;
    w[7] = 1.0*kass_re8*s3 - 1.0*kdiss_re8*s9;
    w[8] = 1.0*kass_re9*s3 - 1.0*kdiss_re9*s7;
    w[9] = 1.0*kass_re10*s4 - 1.0*kdiss_re10*s7;
    w[10] = 1.0*kass_re11*s5 - 1.0*kdiss_re11*s7;
    w[11] = 1.0*kass_re12*s13 - 1.0*kdiss_re12*s14;
    w[12] = 1.0*kass_re13*s7 - 1.0*kdiss_re13*s13;
    w[13] = 1.0*kass_re14*s8 - 1.0*kdiss_re14*s13;
    w[14] = 1.0*kass_re15*s9 - 1.0*kdiss_re15*s13;
    w[15] = 1.0*kass_re16*s10 - 1.0*kdiss_re16*s13;
    w[16] = 1.0*kass_re17*s14 - 1.0*kdiss_re17*s15;
    w[17] = 1.0*kass_re18*s7 - 1.0*kdiss_re18*s15;
    w[18] = 1.0*kass_re19*s14 - 1.0*kdiss_re19*s16;
    w[19] = 1.0*kass_re20*s11 - 1.0*kdiss_re20*s16;
    w[20] = 1.0*kass_re21*s12 - 1.0*kdiss_re21*s16;
    w[21] = 1.0*kass_re22*s17 - 1.0*kdiss_re22*s18;
    w[22] = 1.0*kass_re23*s14 - 1.0*kdiss_re23*s17;
    w[23] = 1.0*kass_re24*s18 - 1.0*kdiss_re24*s19;
    w[24] = 1.0*kass_re25*s18 - 1.0*kdiss_re25*s20;
    w[25] = 1.0*kass_re26*s18 - 1.0*kdiss_re26*s21;
    w[26] = 1.0*kass_re27*s18 - 1.0*kdiss_re27*s22;
    w[27] = 1.0*kass_re28*s18 - 1.0*kdiss_re28*s23;
    w[28] = 1.0*kass_re29*s18 - 1.0*kdiss_re29*s24;
    w[29] = 1.0*kass_re30*s18 - 1.0*kdiss_re30*s25;
    w[30] = 1.0*kass_re31*s18 - 1.0*kdiss_re31*s26;
    w[31] = 1.0*kass_re32*s27 - 1.0*kdiss_re32*s28;
    w[32] = 1.0*kass_re33*s18 - 1.0*kdiss_re33*s27;
    w[33] = 1.0*kass_re34*s15 - 1.0*kdiss_re34*s19;
    w[34] = 1.0*kass_re35*s15 - 1.0*kdiss_re35*s20;
    w[35] = 1.0*kass_re36*s16 - 1.0*kdiss_re36*s26;
    w[36] = 1.0*kass_re37*s28 - 1.0*kdiss_re37*s29;
    w[37] = 1.0*kass_re38*s28 - 1.0*kdiss_re38*s30;
    w[38] = 1.0*kass_re39*s28 - 1.0*kdiss_re39*s31;
    w[39] = 1.0*kass_re40*s28 - 1.0*kdiss_re40*s32;
    w[40] = 1.0*kass_re41*s20 - 1.0*kdiss_re41*s30;
    w[41] = 1.0*kass_re42*s20 - 1.0*kdiss_re42*s31;
    w[42] = 1.0*kass_re43*s20 - 1.0*kdiss_re43*s32;
    w[43] = 1.0*kass_re44*s26 - 1.0*kdiss_re44*s30;
    w[44] = 1.0*kass_re45*s33 - 1.0*kdiss_re45*s34;
    w[45] = 1.0*kass_re46*s35 - 1.0*kdiss_re46*s36;
    w[46] = 1.0*kass_re47*s37 - 1.0*kdiss_re47*s38;
    w[47] = 1.0*kass_re48*s39 - 1.0*kdiss_re48*s40;
    w[48] = 1.0*kass_re49*s41 - 1.0*kdiss_re49*s42;
    w[49] = 1.0*kass_re50*s43 - 1.0*kdiss_re50*s44;
    w[50] = 1.0*kass_re51*s45 - 1.0*kdiss_re51*s46;
    w[51] = 1.0*kass_re52*s47 - 1.0*kdiss_re52*s48;
    w[52] = 1.0*kass_re53*s49 - 1.0*kdiss_re53*s50;
    w[53] = 1.0*kass_re54*s51 - 1.0*kdiss_re54*s52;
    w[54] = 1.0*kass_re55*s29 - 1.0*kdiss_re55*s37;
    w[55] = 1.0*kass_re56*s29 - 1.0*kdiss_re56*s33;
    w[56] = 1.0*kass_re57*s30 - 1.0*kdiss_re57*s35;
    w[57] = 1.0*kass_re58*s30 - 1.0*kdiss_re58*s41;
    w[58] = 1.0*kass_re59*s30 - 1.0*kdiss_re59*s47;
    w[59] = 1.0*kass_re60*s31 - 1.0*kdiss_re60*s33;
    w[60] = 1.0*kass_re61*s31 - 1.0*kdiss_re61*s45;
    w[61] = 1.0*kass_re62*s31 - 1.0*kdiss_re62*s39;
    w[62] = 1.0*kass_re63*s32 - 1.0*kdiss_re63*s47;
    w[63] = 1.0*kass_re64*s32 - 1.0*kdiss_re64*s45;
    w[64] = 1.0*kass_re65*s32 - 1.0*kdiss_re65*s35;
    w[65] = 1.0*kass_re66*s28 - 1.0*kdiss_re66*s56;
    w[66] = 1.0*kass_re67*s28 - 1.0*kdiss_re67*s49;
    w[67] = 1.0*kass_re68*s28 - 1.0*kdiss_re68*s51;
    w[68] = 1.0*kass_re69*s28 - 1.0*kdiss_re69*s53;
    w[69] = 1.0*kass_re70*s28 - 1.0*kdiss_re70*s54;
    w[70] = 1.0*kass_re71*s28 - 1.0*kdiss_re71*s55;
    w[71] = 1.0*kass_re72*s40 - 1.0*kdiss_re72*s57;
    w[72] = 1.0*kass_re73*s53 - 1.0*kdiss_re73*s57;
    w[73] = 1.0*kass_re74*s54 - 1.0*kdiss_re74*s57;
    w[74] = 1.0*kass_re75*s52 - 1.0*kdiss_re75*s57;
    w[75] = 1.0*kass_re76*s50 - 1.0*kdiss_re76*s57;
    w[76] = 1.0*kass_re77*s56 - 1.0*kdiss_re77*s57;
    w[77] = 1.0*kass_re78*s48 - 1.0*kdiss_re78*s57;
    w[78] = 1.0*kass_re79*s30 - 1.0*kdiss_re79*s43;
    w[79] = 1.0*kass_re80*s55 - 1.0*kdiss_re80*s57;
    w[80] = 1.0*kass_re81*s42 - 1.0*kdiss_re81*s57;
    w[81] = 1.0*kass_re82*s44 - 1.0*kdiss_re82*s57;
    w[82] = 1.0*kass_re83*s38 - 1.0*kdiss_re83*s57;
    w[83] = 1.0*kass_re84*s36 - 1.0*kdiss_re84*s57;
    w[84] = 1.0*kass_re85*s34 - 1.0*kdiss_re85*s57;
    w[85] = 1.0*kass_re86*s46 - 1.0*kdiss_re86*s57;
}