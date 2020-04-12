#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_orfao1(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 5;
    rowvals[2] = 3;
    rowvals[3] = 6;
    rowvals[4] = 7;
    rowvals[5] = 2;
    rowvals[6] = 2;
    rowvals[7] = 4;
    rowvals[8] = 0;
    rowvals[9] = 3;
    rowvals[10] = 2;
    rowvals[11] = 0;
    rowvals[12] = 1;
    rowvals[13] = 2;
    rowvals[14] = 5;
}