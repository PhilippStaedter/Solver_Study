#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_borghans3(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 6;
    rowvals[2] = 2;
    rowvals[3] = 4;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 2;
    rowvals[7] = 5;
}