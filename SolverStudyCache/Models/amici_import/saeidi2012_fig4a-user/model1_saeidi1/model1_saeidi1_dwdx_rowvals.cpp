#include "sundials/sundials_types.h"

void dwdx_rowvals_model1_saeidi1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 5;
    rowvals[2] = 1;
    rowvals[3] = 5;
    rowvals[4] = 3;
    rowvals[5] = 4;
    rowvals[6] = 0;
    rowvals[7] = 2;
    rowvals[8] = 4;
    rowvals[9] = 6;
    rowvals[10] = 4;
    rowvals[11] = 5;
    rowvals[12] = 1;
}