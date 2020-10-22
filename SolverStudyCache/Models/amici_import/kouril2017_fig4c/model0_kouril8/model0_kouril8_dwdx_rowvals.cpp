#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_kouril8(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 0;
    rowvals[4] = 2;
    rowvals[5] = 3;
    rowvals[6] = 0;
    rowvals[7] = 3;
    rowvals[8] = 5;
    rowvals[9] = 4;
    rowvals[10] = 3;
    rowvals[11] = 4;
    rowvals[12] = 3;
    rowvals[13] = 1;
    rowvals[14] = 3;
}