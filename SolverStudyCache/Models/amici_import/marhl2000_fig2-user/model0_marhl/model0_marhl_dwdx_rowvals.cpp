#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_marhl(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 3;
    rowvals[2] = 5;
    rowvals[3] = 1;
    rowvals[4] = 0;
    rowvals[5] = 2;
    rowvals[6] = 3;
    rowvals[7] = 4;
    rowvals[8] = 5;
    rowvals[9] = 6;
    rowvals[10] = 2;
}