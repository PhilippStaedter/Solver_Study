cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(fung1_Fig3A_Vgly_0_5)

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

set(SRC_LIST_LIB fung1_Fig3A_Vgly_0_5_J.cpp 
fung1_Fig3A_Vgly_0_5_JB.cpp 
fung1_Fig3A_Vgly_0_5_JDiag.cpp 
fung1_Fig3A_Vgly_0_5_JSparse.cpp 
fung1_Fig3A_Vgly_0_5_JSparseB.cpp 
fung1_Fig3A_Vgly_0_5_Jy.cpp 
fung1_Fig3A_Vgly_0_5_dJydsigmay.cpp 
fung1_Fig3A_Vgly_0_5_dJydy.cpp 
fung1_Fig3A_Vgly_0_5_dwdp.cpp 
fung1_Fig3A_Vgly_0_5_dwdx.cpp 
fung1_Fig3A_Vgly_0_5_dxdotdw.cpp 
fung1_Fig3A_Vgly_0_5_dxdotdp.cpp 
fung1_Fig3A_Vgly_0_5_dydx.cpp 
fung1_Fig3A_Vgly_0_5_dydp.cpp 
fung1_Fig3A_Vgly_0_5_dsigmaydp.cpp 
fung1_Fig3A_Vgly_0_5_sigmay.cpp 
fung1_Fig3A_Vgly_0_5_w.cpp 
fung1_Fig3A_Vgly_0_5_x0.cpp 
fung1_Fig3A_Vgly_0_5_x0_fixedParameters.cpp 
fung1_Fig3A_Vgly_0_5_sx0.cpp 
fung1_Fig3A_Vgly_0_5_sx0_fixedParameters.cpp 
fung1_Fig3A_Vgly_0_5_xdot.cpp 
fung1_Fig3A_Vgly_0_5_y.cpp 
fung1_Fig3A_Vgly_0_5_x_rdata.cpp 
fung1_Fig3A_Vgly_0_5_total_cl.cpp 
fung1_Fig3A_Vgly_0_5_x_solver.cpp 
fung1_Fig3A_Vgly_0_5_JSparse_colptrs.cpp
fung1_Fig3A_Vgly_0_5_JSparse_rowvals.cpp 
fung1_Fig3A_Vgly_0_5_JSparseB_colptrs.cpp
fung1_Fig3A_Vgly_0_5_JSparseB_rowvals.cpp 
fung1_Fig3A_Vgly_0_5_dJydy_colptrs.cpp
fung1_Fig3A_Vgly_0_5_dJydy_rowvals.cpp 
fung1_Fig3A_Vgly_0_5_dwdx_colptrs.cpp
fung1_Fig3A_Vgly_0_5_dwdx_rowvals.cpp 
fung1_Fig3A_Vgly_0_5_dxdotdw_colptrs.cpp
fung1_Fig3A_Vgly_0_5_dxdotdw_rowvals.cpp 
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

