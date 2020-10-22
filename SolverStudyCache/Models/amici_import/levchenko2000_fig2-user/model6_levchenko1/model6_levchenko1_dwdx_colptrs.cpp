#include "sundials/sundials_types.h"

void dwdx_colptrs_model6_levchenko1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 5;
    colptrs[4] = 7;
    colptrs[5] = 9;
    colptrs[6] = 11;
    colptrs[7] = 12;
    colptrs[8] = 14;
    colptrs[9] = 15;
    colptrs[10] = 17;
    colptrs[11] = 19;
    colptrs[12] = 21;
    colptrs[13] = 23;
    colptrs[14] = 25;
    colptrs[15] = 28;
    colptrs[16] = 30;
    colptrs[17] = 31;
    colptrs[18] = 32;
    colptrs[19] = 33;
    colptrs[20] = 35;
    colptrs[21] = 38;
    colptrs[22] = 40;
}