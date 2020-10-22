#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_piedrafita1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 4;
    colptrs[3] = 9;
    colptrs[4] = 12;
    colptrs[5] = 14;
    colptrs[6] = 16;
    colptrs[7] = 20;
    colptrs[8] = 22;
    colptrs[9] = 24;
    colptrs[10] = 25;
    colptrs[11] = 27;
}