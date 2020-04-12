#include "sundials/sundials_types.h"

void dwdx_rowvals_model3_chance3(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 1;
    rowvals[3] = 0;
}