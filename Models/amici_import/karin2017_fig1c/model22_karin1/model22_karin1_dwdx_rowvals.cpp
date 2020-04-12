#include "sundials/sundials_types.h"

void dwdx_rowvals_model22_karin1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 2;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 2;
}