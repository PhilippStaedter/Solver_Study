#include "sundials/sundials_types.h"

void dwdx_rowvals_model2_alsheihk1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 9;
    rowvals[2] = 2;
    rowvals[3] = 4;
    rowvals[4] = 5;
    rowvals[5] = 6;
    rowvals[6] = 2;
    rowvals[7] = 7;
    rowvals[8] = 8;
    rowvals[9] = 2;
    rowvals[10] = 3;
}