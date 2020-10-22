cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(Liu2011)

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

set(SRC_LIST_LIB Liu2011_J.cpp 
Liu2011_JB.cpp 
Liu2011_JDiag.cpp 
Liu2011_JSparse.cpp 
Liu2011_JSparseB.cpp 
Liu2011_Jy.cpp 
Liu2011_dJydsigmay.cpp 
Liu2011_dJydy.cpp 
Liu2011_dwdp.cpp 
Liu2011_dwdx.cpp 
Liu2011_dxdotdw.cpp 
Liu2011_dxdotdp.cpp 
Liu2011_dydx.cpp 
Liu2011_dydp.cpp 
Liu2011_dsigmaydp.cpp 
Liu2011_sigmay.cpp 
Liu2011_w.cpp 
Liu2011_x0.cpp 
Liu2011_x0_fixedParameters.cpp 
Liu2011_sx0.cpp 
Liu2011_sx0_fixedParameters.cpp 
Liu2011_xdot.cpp 
Liu2011_y.cpp 
Liu2011_x_rdata.cpp 
Liu2011_total_cl.cpp 
Liu2011_x_solver.cpp 
Liu2011_JSparse_colptrs.cpp
Liu2011_JSparse_rowvals.cpp 
Liu2011_JSparseB_colptrs.cpp
Liu2011_JSparseB_rowvals.cpp 
Liu2011_dJydy_colptrs.cpp
Liu2011_dJydy_rowvals.cpp 
Liu2011_dwdx_colptrs.cpp
Liu2011_dwdx_rowvals.cpp 
Liu2011_dxdotdw_colptrs.cpp
Liu2011_dxdotdw_rowvals.cpp 
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

