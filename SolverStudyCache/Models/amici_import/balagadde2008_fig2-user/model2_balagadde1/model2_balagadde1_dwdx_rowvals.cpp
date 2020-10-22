#include "sundials/sundials_types.h"

void dwdx_rowvals_model2_balagadde1(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 5;
    rowvals[2] = 1;
    rowvals[3] = 7;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 2;
    rowvals[7] = 4;
    rowvals[8] = 0;
    rowvals[9] = 2;
    rowvals[10] = 3;
    rowvals[11] = 6;
    rowvals[12] = 1;
    rowvals[13] = 6;
}