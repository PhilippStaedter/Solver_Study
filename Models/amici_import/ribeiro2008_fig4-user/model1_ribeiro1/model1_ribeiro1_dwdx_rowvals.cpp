#include "sundials/sundials_types.h"

void dwdx_rowvals_model1_ribeiro1(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 4;
    rowvals[2] = 1;
}