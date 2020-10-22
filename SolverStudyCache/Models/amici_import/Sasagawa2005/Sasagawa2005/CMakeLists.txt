cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(Sasagawa2005)

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

set(SRC_LIST_LIB Sasagawa2005_J.cpp 
Sasagawa2005_JB.cpp 
Sasagawa2005_JDiag.cpp 
Sasagawa2005_JSparse.cpp 
Sasagawa2005_JSparseB.cpp 
Sasagawa2005_Jy.cpp 
Sasagawa2005_dJydsigmay.cpp 
Sasagawa2005_dJydy.cpp 
Sasagawa2005_dwdp.cpp 
Sasagawa2005_dwdx.cpp 
Sasagawa2005_dxdotdw.cpp 
Sasagawa2005_dxdotdp.cpp 
Sasagawa2005_dydx.cpp 
Sasagawa2005_dydp.cpp 
Sasagawa2005_dsigmaydp.cpp 
Sasagawa2005_sigmay.cpp 
Sasagawa2005_w.cpp 
Sasagawa2005_x0.cpp 
Sasagawa2005_x0_fixedParameters.cpp 
Sasagawa2005_sx0.cpp 
Sasagawa2005_sx0_fixedParameters.cpp 
Sasagawa2005_xdot.cpp 
Sasagawa2005_y.cpp 
Sasagawa2005_x_rdata.cpp 
Sasagawa2005_total_cl.cpp 
Sasagawa2005_x_solver.cpp 
Sasagawa2005_JSparse_colptrs.cpp
Sasagawa2005_JSparse_rowvals.cpp 
Sasagawa2005_JSparseB_colptrs.cpp
Sasagawa2005_JSparseB_rowvals.cpp 
Sasagawa2005_dJydy_colptrs.cpp
Sasagawa2005_dJydy_rowvals.cpp 
Sasagawa2005_dwdx_colptrs.cpp
Sasagawa2005_dwdx_rowvals.cpp 
Sasagawa2005_dxdotdw_colptrs.cpp
Sasagawa2005_dxdotdw_rowvals.cpp 
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

