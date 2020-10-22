#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_panteleev1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 10;
    colptrs[5] = 13;
    colptrs[6] = 15;
    colptrs[7] = 17;
    colptrs[8] = 20;
    colptrs[9] = 22;
}