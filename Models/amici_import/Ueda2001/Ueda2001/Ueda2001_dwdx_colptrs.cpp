#include "sundials/sundials_types.h"

void dwdx_colptrs_Ueda2001(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 4;
    colptrs[3] = 10;
    colptrs[4] = 13;
    colptrs[5] = 16;
    colptrs[6] = 19;
    colptrs[7] = 22;
    colptrs[8] = 26;
    colptrs[9] = 32;
    colptrs[10] = 35;
    colptrs[11] = 38;
    colptrs[12] = 39;
    colptrs[13] = 40;
}