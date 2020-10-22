#include "amici/model.h"
#include "wrapfunctions.h"

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(new Model_fung1_Fig3A_Vgly_0_5());
}
