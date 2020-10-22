#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_kholodenko2(sunindextype *rowvals){
    rowvals[0] = 6;
    rowvals[1] = 7;
    rowvals[2] = 9;
    rowvals[3] = 0;
    rowvals[4] = 8;
    rowvals[5] = 2;
    rowvals[6] = 0;
    rowvals[7] = 1;
    rowvals[8] = 2;
    rowvals[9] = 3;
    rowvals[10] = 3;
    rowvals[11] = 5;
    rowvals[12] = 4;
    rowvals[13] = 6;
    rowvals[14] = 7;
}