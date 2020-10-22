#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_miao2008(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 1;
    rowvals[5] = 4;
    rowvals[6] = 5;
    rowvals[7] = 7;
    rowvals[8] = 3;
    rowvals[9] = 8;
    rowvals[10] = 2;
    rowvals[11] = 5;
    rowvals[12] = 6;
    rowvals[13] = 7;
}