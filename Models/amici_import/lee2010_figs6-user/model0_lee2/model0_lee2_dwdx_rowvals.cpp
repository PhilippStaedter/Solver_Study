#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_lee2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 4;
    rowvals[3] = 6;
    rowvals[4] = 4;
    rowvals[5] = 5;
    rowvals[6] = 1;
    rowvals[7] = 2;
    rowvals[8] = 0;
    rowvals[9] = 3;
    rowvals[10] = 6;
    rowvals[11] = 7;
    rowvals[12] = 4;
    rowvals[13] = 0;
    rowvals[14] = 6;
    rowvals[15] = 1;
}