#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_zhao1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 5;
    rowvals[3] = 9;
    rowvals[4] = 2;
    rowvals[5] = 5;
    rowvals[6] = 7;
    rowvals[7] = 8;
    rowvals[8] = 2;
    rowvals[9] = 5;
    rowvals[10] = 6;
    rowvals[11] = 2;
    rowvals[12] = 3;
    rowvals[13] = 5;
}