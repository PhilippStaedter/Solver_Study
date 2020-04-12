cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(model1_alexander1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

find_package(Amici 0.10.7 HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)
message(STATUS "Found AMICI ${Amici_DIR}")

set(MODEL_DIR ${CMAKE_CURRENT_LIST_DIR})

set(SRC_LIST_LIB model1_alexander1_J.cpp 
model1_alexander1_JB.cpp 
model1_alexander1_JDiag.cpp 
model1_alexander1_JSparse.cpp 
model1_alexander1_JSparseB.cpp 
model1_alexander1_Jy.cpp 
model1_alexander1_dJydsigmay.cpp 
model1_alexander1_dJydy.cpp 
model1_alexander1_dwdp.cpp 
model1_alexander1_dwdx.cpp 
model1_alexander1_dxdotdw.cpp 
model1_alexander1_dxdotdp.cpp 
model1_alexander1_dydx.cpp 
model1_alexander1_dydp.cpp 
model1_alexander1_dsigmaydp.cpp 
model1_alexander1_sigmay.cpp 
model1_alexander1_w.cpp 
model1_alexander1_x0.cpp 
model1_alexander1_x0_fixedParameters.cpp 
model1_alexander1_sx0.cpp 
model1_alexander1_sx0_fixedParameters.cpp 
model1_alexander1_xdot.cpp 
model1_alexander1_y.cpp 
model1_alexander1_x_rdata.cpp 
model1_alexander1_total_cl.cpp 
model1_alexander1_x_solver.cpp 
model1_alexander1_JSparse_colptrs.cpp
model1_alexander1_JSparse_rowvals.cpp 
model1_alexander1_JSparseB_colptrs.cpp
model1_alexander1_JSparseB_rowvals.cpp 
model1_alexander1_dJydy_colptrs.cpp
model1_alexander1_dJydy_rowvals.cpp 
model1_alexander1_dwdx_colptrs.cpp
model1_alexander1_dwdx_rowvals.cpp 
model1_alexander1_dxdotdw_colptrs.cpp
model1_alexander1_dxdotdw_rowvals.cpp 
${MODEL_DIR}/wrapfunctions.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
add_library(model ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(${PROJECT_NAME}
    PUBLIC Upstream::amici
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

## SWIG
option(ENABLE_SWIG "Build swig/python library?" ON)
if(ENABLE_SWIG)
    if(NOT(${CMAKE_VERSION} VERSION_LESS 3.8))
        add_subdirectory(swig)
    else()
        message(WARNING "Unable to build SWIG interface, upgrade CMake to >=3.8.")
    endif()
endif()


# <Export cmake configuration>
include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
export(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE Upstream::
    )
# </Export cmake configuration>

