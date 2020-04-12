#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_karin5(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 3;
    rowvals[3] = 0;
    rowvals[4] = 1;
    rowvals[5] = 2;
    rowvals[6] = 3;
    rowvals[7] = 0;
    rowvals[8] = 1;
    rowvals[9] = 1;
    rowvals[10] = 3;
    rowvals[11] = 2;
    rowvals[12] = 3;
    rowvals[13] = 4;
}