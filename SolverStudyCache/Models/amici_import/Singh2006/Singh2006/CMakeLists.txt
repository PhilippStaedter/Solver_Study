cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(Singh2006)

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

set(SRC_LIST_LIB Singh2006_J.cpp 
Singh2006_JB.cpp 
Singh2006_JDiag.cpp 
Singh2006_JSparse.cpp 
Singh2006_JSparseB.cpp 
Singh2006_Jy.cpp 
Singh2006_dJydsigmay.cpp 
Singh2006_dJydy.cpp 
Singh2006_dwdp.cpp 
Singh2006_dwdx.cpp 
Singh2006_dxdotdw.cpp 
Singh2006_dxdotdp.cpp 
Singh2006_dydx.cpp 
Singh2006_dydp.cpp 
Singh2006_dsigmaydp.cpp 
Singh2006_sigmay.cpp 
Singh2006_w.cpp 
Singh2006_x0.cpp 
Singh2006_x0_fixedParameters.cpp 
Singh2006_sx0.cpp 
Singh2006_sx0_fixedParameters.cpp 
Singh2006_xdot.cpp 
Singh2006_y.cpp 
Singh2006_x_rdata.cpp 
Singh2006_total_cl.cpp 
Singh2006_x_solver.cpp 
Singh2006_JSparse_colptrs.cpp
Singh2006_JSparse_rowvals.cpp 
Singh2006_JSparseB_colptrs.cpp
Singh2006_JSparseB_rowvals.cpp 
Singh2006_dJydy_colptrs.cpp
Singh2006_dJydy_rowvals.cpp 
Singh2006_dwdx_colptrs.cpp
Singh2006_dwdx_rowvals.cpp 
Singh2006_dxdotdw_colptrs.cpp
Singh2006_dxdotdw_rowvals.cpp 
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

