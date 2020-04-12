#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_panteleev2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 6;
    colptrs[3] = 8;
    colptrs[4] = 10;
    colptrs[5] = 13;
    colptrs[6] = 14;
    colptrs[7] = 16;
    colptrs[8] = 19;
    colptrs[9] = 21;
}