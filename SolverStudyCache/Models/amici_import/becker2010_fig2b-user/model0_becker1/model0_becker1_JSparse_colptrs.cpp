#include "sundials/sundials_types.h"

void JSparse_colptrs_model0_becker1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 10;
    colptrs[4] = 15;
    colptrs[5] = 15;
    colptrs[6] = 15;
}