#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_kim5(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 7;
    rowvals[2] = 1;
    rowvals[3] = 6;
    rowvals[4] = 2;
    rowvals[5] = 4;
    rowvals[6] = 3;
    rowvals[7] = 5;
    rowvals[8] = 8;
}