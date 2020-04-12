#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_panteleev3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 7;
    colptrs[5] = 8;
    colptrs[6] = 10;
    colptrs[7] = 12;
    colptrs[8] = 13;
}