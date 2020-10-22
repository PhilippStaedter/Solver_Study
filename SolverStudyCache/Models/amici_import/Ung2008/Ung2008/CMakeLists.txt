cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(Ung2008)

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

set(SRC_LIST_LIB Ung2008_J.cpp 
Ung2008_JB.cpp 
Ung2008_JDiag.cpp 
Ung2008_JSparse.cpp 
Ung2008_JSparseB.cpp 
Ung2008_Jy.cpp 
Ung2008_dJydsigmay.cpp 
Ung2008_dJydy.cpp 
Ung2008_dwdp.cpp 
Ung2008_dwdx.cpp 
Ung2008_dxdotdw.cpp 
Ung2008_dxdotdp.cpp 
Ung2008_dydx.cpp 
Ung2008_dydp.cpp 
Ung2008_dsigmaydp.cpp 
Ung2008_sigmay.cpp 
Ung2008_w.cpp 
Ung2008_x0.cpp 
Ung2008_x0_fixedParameters.cpp 
Ung2008_sx0.cpp 
Ung2008_sx0_fixedParameters.cpp 
Ung2008_xdot.cpp 
Ung2008_y.cpp 
Ung2008_x_rdata.cpp 
Ung2008_total_cl.cpp 
Ung2008_x_solver.cpp 
Ung2008_JSparse_colptrs.cpp
Ung2008_JSparse_rowvals.cpp 
Ung2008_JSparseB_colptrs.cpp
Ung2008_JSparseB_rowvals.cpp 
Ung2008_dJydy_colptrs.cpp
Ung2008_dJydy_rowvals.cpp 
Ung2008_dwdx_colptrs.cpp
Ung2008_dwdx_rowvals.cpp 
Ung2008_dxdotdw_colptrs.cpp
Ung2008_dxdotdw_rowvals.cpp 
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

