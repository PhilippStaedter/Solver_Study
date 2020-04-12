#include "sundials/sundials_types.h"

void dwdx_colptrs_Leber2015(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 8;
    colptrs[2] = 12;
    colptrs[3] = 16;
    colptrs[4] = 17;
    colptrs[5] = 19;
    colptrs[6] = 22;
    colptrs[7] = 27;
    colptrs[8] = 31;
    colptrs[9] = 33;
    colptrs[10] = 33;
    colptrs[11] = 38;
    colptrs[12] = 42;
    colptrs[13] = 43;
    colptrs[14] = 44;
    colptrs[15] = 46;
    colptrs[16] = 52;
    colptrs[17] = 53;
    colptrs[18] = 57;
    colptrs[19] = 60;
    colptrs[20] = 61;
    colptrs[21] = 61;
    colptrs[22] = 62;
    colptrs[23] = 63;
}