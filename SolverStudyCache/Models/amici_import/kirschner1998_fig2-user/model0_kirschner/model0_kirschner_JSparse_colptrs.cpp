#include "sundials/sundials_types.h"

void JSparse_colptrs_model0_kirschner(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
}