#include "sundials/sundials_types.h"

void JSparse_colptrs_lou1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 5;
    colptrs[2] = 10;
    colptrs[3] = 16;
    colptrs[4] = 22;
    colptrs[5] = 28;
    colptrs[6] = 34;
}