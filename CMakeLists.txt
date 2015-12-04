project(FractureDynamics)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
set(IFEM_PRECONFIGURED ${IFEM_CONFIGURED})
if (NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
endif(NOT IFEM_CONFIGURED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR}
                    ../Elasticity ../Elasticity/Linear)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif(NOT WIN32)

file(GLOB LinEl_SRCS ../Elasticity/*.C)

add_library(CommonCH CahnHilliard.C FractureElasticity.C ${LinEl_SRCS})

add_executable(CahnHilliard main_CahnHilliard.C)
add_executable(FractureDynamics main_FractureDynamics.C)

target_link_libraries(CahnHilliard CommonCH ${IFEM_LIBRARIES})
target_link_libraries(FractureDynamics CommonCH ${IFEM_LIBRARIES})

# Installation
install(TARGETS CahnHilliard FractureDynamics DESTINATION bin)

# Testing
enable_testing()

# Generate regtest script with correct paths
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

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
                  CommonCH ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
