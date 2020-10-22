#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model2_essunger3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 10;
    colptrs[6] = 11;
    colptrs[7] = 13;
    colptrs[8] = 15;
    colptrs[9] = 17;
    colptrs[10] = 19;
    colptrs[11] = 21;
    colptrs[12] = 23;
    colptrs[13] = 25;
    colptrs[14] = 27;
    colptrs[15] = 28;
    colptrs[16] = 29;
    colptrs[17] = 31;
}