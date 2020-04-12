#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_bachar1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 3;
    rowvals[3] = 4;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 5;
    rowvals[7] = 6;
    rowvals[8] = 7;
    rowvals[9] = 0;
    rowvals[10] = 1;
    rowvals[11] = 2;
}