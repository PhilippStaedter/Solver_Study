cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(model24_karin3)

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

set(SRC_LIST_LIB model24_karin3_J.cpp 
model24_karin3_JB.cpp 
model24_karin3_JDiag.cpp 
model24_karin3_JSparse.cpp 
model24_karin3_JSparseB.cpp 
model24_karin3_Jy.cpp 
model24_karin3_dJydsigmay.cpp 
model24_karin3_dJydy.cpp 
model24_karin3_dwdp.cpp 
model24_karin3_dwdx.cpp 
model24_karin3_dxdotdw.cpp 
model24_karin3_dxdotdp.cpp 
model24_karin3_dydx.cpp 
model24_karin3_dydp.cpp 
model24_karin3_dsigmaydp.cpp 
model24_karin3_sigmay.cpp 
model24_karin3_w.cpp 
model24_karin3_x0.cpp 
model24_karin3_x0_fixedParameters.cpp 
model24_karin3_sx0.cpp 
model24_karin3_sx0_fixedParameters.cpp 
model24_karin3_xdot.cpp 
model24_karin3_y.cpp 
model24_karin3_x_rdata.cpp 
model24_karin3_total_cl.cpp 
model24_karin3_x_solver.cpp 
model24_karin3_JSparse_colptrs.cpp
model24_karin3_JSparse_rowvals.cpp 
model24_karin3_JSparseB_colptrs.cpp
model24_karin3_JSparseB_rowvals.cpp 
model24_karin3_dJydy_colptrs.cpp
model24_karin3_dJydy_rowvals.cpp 
model24_karin3_dwdx_colptrs.cpp
model24_karin3_dwdx_rowvals.cpp 
model24_karin3_dxdotdw_colptrs.cpp
model24_karin3_dxdotdw_rowvals.cpp 
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

