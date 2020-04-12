cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(Sivakumar2011c)

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

set(SRC_LIST_LIB Sivakumar2011c_J.cpp 
Sivakumar2011c_JB.cpp 
Sivakumar2011c_JDiag.cpp 
Sivakumar2011c_JSparse.cpp 
Sivakumar2011c_JSparseB.cpp 
Sivakumar2011c_Jy.cpp 
Sivakumar2011c_dJydsigmay.cpp 
Sivakumar2011c_dJydy.cpp 
Sivakumar2011c_dwdp.cpp 
Sivakumar2011c_dwdx.cpp 
Sivakumar2011c_dxdotdw.cpp 
Sivakumar2011c_dxdotdp.cpp 
Sivakumar2011c_dydx.cpp 
Sivakumar2011c_dydp.cpp 
Sivakumar2011c_dsigmaydp.cpp 
Sivakumar2011c_sigmay.cpp 
Sivakumar2011c_w.cpp 
Sivakumar2011c_x0.cpp 
Sivakumar2011c_x0_fixedParameters.cpp 
Sivakumar2011c_sx0.cpp 
Sivakumar2011c_sx0_fixedParameters.cpp 
Sivakumar2011c_xdot.cpp 
Sivakumar2011c_y.cpp 
Sivakumar2011c_x_rdata.cpp 
Sivakumar2011c_total_cl.cpp 
Sivakumar2011c_x_solver.cpp 
Sivakumar2011c_JSparse_colptrs.cpp
Sivakumar2011c_JSparse_rowvals.cpp 
Sivakumar2011c_JSparseB_colptrs.cpp
Sivakumar2011c_JSparseB_rowvals.cpp 
Sivakumar2011c_dJydy_colptrs.cpp
Sivakumar2011c_dJydy_rowvals.cpp 
Sivakumar2011c_dwdx_colptrs.cpp
Sivakumar2011c_dwdx_rowvals.cpp 
Sivakumar2011c_dxdotdw_colptrs.cpp
Sivakumar2011c_dxdotdw_rowvals.cpp 
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

