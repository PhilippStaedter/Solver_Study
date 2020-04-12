#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_mellor1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 3;
    colptrs[3] = 3;
    colptrs[4] = 3;
    colptrs[5] = 3;
    colptrs[6] = 3;
    colptrs[7] = 3;
    colptrs[8] = 4;
    colptrs[9] = 5;
    colptrs[10] = 5;
}