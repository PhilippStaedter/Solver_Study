#include "sundials/sundials_types.h"

void dwdx_rowvals_model4_mellor1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 4;
}