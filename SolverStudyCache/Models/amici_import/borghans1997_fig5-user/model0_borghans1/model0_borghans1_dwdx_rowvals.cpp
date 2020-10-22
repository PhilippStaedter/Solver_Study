#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_borghans1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 4;
    rowvals[2] = 5;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 0;
    rowvals[6] = 1;
    rowvals[7] = 3;
    rowvals[8] = 4;
}