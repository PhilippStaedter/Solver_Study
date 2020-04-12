#include "sundials/sundials_types.h"

void dwdx_rowvals_kolodkin2(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 0;
    rowvals[4] = 2;
    rowvals[5] = 2;
    rowvals[6] = 0;
    rowvals[7] = 0;
}