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

if(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif(NOT WIN32)

if(EXISTS ${PROJECT_SOURCE_DIR}/../IFEM-PoroElasticity)
  if(NOT TARGET PoroElastic)
    add_subdirectory(../IFEM-PoroElasticity/PoroElastic PoroElasticity)
  endif()
  include_directories(../IFEM-PoroElasticity/PoroElastic)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/../PoroElasticity)
  if(NOT TARGET PoroElastic)
    add_subdirectory(../PoroElasticity/PoroElastic PoroElasticity)
  endif()
  include_directories(../PoroElasticity/PoroElastic)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/../IFEM-Elasticity)
  if(NOT TARGET Elasticity)
    add_subdirectory(../IFEM-Elasticity Elasticity)
  endif()
  include_directories(../IFEM-Elasticity)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/../Elasticity)
  if(NOT TARGET Elasticity)
    add_subdirectory(../Elasticity Elasticity)
  endif()
  include_directories(../Elasticity)
endif()

file(GLOB FracEl_SRCS [A-Z]*.C)

add_library(CommonFrac STATIC ${FracEl_SRCS})
set(Common_LIBRARIES CommonFrac Elasticity ${IFEM_LIBRARIES})
set(Fracture_LIBRARIES CommonFrac PoroElastic Elasticity ${IFEM_LIBRARIES})

add_executable(CahnHilliard main_CahnHilliard.C)
add_executable(FractureDynamics main_FractureDynamics.C)

target_link_libraries(CahnHilliard ${Common_LIBRARIES})
target_link_libraries(FractureDynamics ${Fracture_LIBRARIES})

# Installation
install(TARGETS CahnHilliard FractureDynamics DESTINATION bin)

# Testing
enable_testing()
include(IFEMTesting)

# Generate regtest script with correct paths
configure_file(${IFEM_REGTEST_SCRIPT} regtest.sh)

# Regression tests
if(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
else()
  file(GLOB CH_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                         ${PROJECT_SOURCE_DIR}/Test/*_CH.reg)
  file(GLOB FD_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                         ${PROJECT_SOURCE_DIR}/Test/*_FD.reg)
  foreach(TESTFILE ${CH_TESTFILES})
    ifem_add_test(${TESTFILE} CahnHilliard)
  endforeach()
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
