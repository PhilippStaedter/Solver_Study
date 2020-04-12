#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_borghans2(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 5;
    rowvals[2] = 6;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 0;
    rowvals[6] = 1;
    rowvals[7] = 3;
    rowvals[8] = 5;
}