#include "amici/model.h"
#include "wrapfunctions.h"

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(new Model_model0_kholodenko1());
}
