project(FractureDynamics)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")

if(EXISTS ${PROJECT_SOURCE_DIR}/../IFEM-Elasticity)
  set(ELASTICITY_DIR ../IFEM-Elasticity)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/../Elasticity)
  set(ELASTICITY_DIR ../Elasticity)
else()
  message(FATAL_ERROR "Need IFEM-Elasticity in a sibling directory.")
endif()

set(EXTRA_DOXY_PATHS "${EXTRA_DOXY_PATHS} \\
                      ${PROJECT_SOURCE_DIR}/${ELASTICITY_DIR}")
if(NOT TARGET Elasticity)
  add_subdirectory(${ELASTICITY_DIR} Elasticity)
endif()
include_directories(${ELASTICITY_DIR})

set(POROELASTIC_DIR ${PROJECT_SOURCE_DIR}/../IFEM-PoroElasticity/PoroElastic)
if(NOT EXISTS ${POROELASTIC_DIR})
  set(POROELASTIC_DIR ${PROJECT_SOURCE_DIR}/../PoroElasticity/PoroElastic)
endif()

if(EXISTS ${POROELASTIC_DIR})
  set(EXTRA_DOXY_PATHS "${EXTRA_DOXY_PATHS} ${POROELASTIC_DIR}")
  if(NOT TARGET PoroElastic)
    add_subdirectory(${POROELASTIC_DIR} PoroElasticity)
  endif()
  include_directories(${POROELASTIC_DIR})
endif()

file(GLOB FracEl_SRCS [A-Z]*.C)

add_library(CommonFrac STATIC ${FracEl_SRCS})
set(Common_LIBRARIES CommonFrac Elasticity ${IFEM_LIBRARIES})
if(EXISTS ${POROELASTIC_DIR})
  set(Fracture_LIBRARIES CommonFrac PoroElastic Elasticity ${IFEM_LIBRARIES})
else()
  set(Fracture_LIBRARIES ${Common_LIBRARIES})
endif()

add_executable(CahnHilliard main_CahnHilliard.C)
add_executable(FractureDynamics main_FractureDynamics.C)

target_link_libraries(CahnHilliard ${Common_LIBRARIES})
target_link_libraries(FractureDynamics ${Fracture_LIBRARIES})

list(APPEND CHECK_SOURCES main_CahnHilliard.C main_FractureDynamics.C ${FracEl_SRCS})

# Installation
install(TARGETS CahnHilliard FractureDynamics DESTINATION bin)

# For generating the doxy
add_doc_target(OpenFrac FractureDynamics)

# Testing
enable_testing()
include(IFEMTesting)

# Regression tests
if(MPI_FOUND)
  # Add parallel tests here
else()
  file(GLOB CH_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                         ${PROJECT_SOURCE_DIR}/Test/*_CH.reg)
  file(GLOB FD_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                         ${PROJECT_SOURCE_DIR}/Test/*_FD.reg)
  foreach(TESTFILE ${CH_TESTFILES})
    ifem_add_test(${TESTFILE} CahnHilliard)
  endforeach()
  # Clang 3.8 optimizer bug trigger with this test
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 3.8)
    list(REMOVE_ITEM FD_TESTFILES Short10x20-adap_FD.reg)
  endif()
  if(LRSpline_FOUND)
    ifem_add_vtf_test(Short10x20-adap_FD.vreg FractureDynamics)
  else()
    list(REMOVE_ITEM FD_TESTFILES Short10x20-adap_FD.reg)
  endif()
  foreach(TESTFILE ${FD_TESTFILES})
    ifem_add_test(${TESTFILE} FractureDynamics)
  endforeach()
endif()
list(APPEND TEST_APPS CahnHilliard FractureDynamics)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  FractureDynamics
                  ${Fracture_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
