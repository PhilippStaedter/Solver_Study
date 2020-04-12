#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_sarma1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 7;
    colptrs[4] = 9;
    colptrs[5] = 15;
    colptrs[6] = 17;
    colptrs[7] = 18;
    colptrs[8] = 20;
    colptrs[9] = 22;
    colptrs[10] = 23;
    colptrs[11] = 25;
    colptrs[12] = 29;
    colptrs[13] = 34;
    colptrs[14] = 36;
    colptrs[15] = 38;
    colptrs[16] = 39;
    colptrs[17] = 41;
    colptrs[18] = 43;
    colptrs[19] = 44;
    colptrs[20] = 45;
    colptrs[21] = 47;
    colptrs[22] = 49;
    colptrs[23] = 51;
    colptrs[24] = 52;
    colptrs[25] = 55;
    colptrs[26] = 58;
    colptrs[27] = 60;
}