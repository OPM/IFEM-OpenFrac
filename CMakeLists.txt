cmake_minimum_required(VERSION 3.22)

project(FractureDynamics LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

ifem_add_dependency_dir(APP Elasticity TARGET Elasticity)
ifem_add_dependency_dir(APP PoroElasticity DIR PoroElastic TARGET PoroElastic OPTIONAL)

ifem_add_library(
  NAME
    CommonFrac
  SOURCES
    CahnHilliard.C
    FractureElasticity.C
    FractureElasticityVoigt.C
    SIMExplPhaseField.C
    SIMDynElasticity.C
    SIMFractureDynamics.C
    SIMFractureQstatic.C
    SIMPhaseField.C
  HEADERS
    CahnHilliard.h
    FractureElasticity.h
    FractureElasticityVoigt.h
    SIMDynElasticity.h
    SIMExplPhaseField.h
    SIMFractureDynamics.h
    SIMFractureQstatic.h
    SIMPhaseField.h
  LIBRARIES
    IFEMAppCommon
    Elasticity
)

if(PoroElastic_DIRECTORY)
  target_compile_definitions(CommonFrac PUBLIC IFEM_HAS_POROELASTIC=1)
  target_sources(CommonFrac PRIVATE PoroFracture.C PoroFracture.h)
  target_link_libraries(CommonFrac PUBLIC PoroElastic)
endif()

ifem_add_application(
  NAME
    CahnHilliard
  SOURCES
    main_CahnHilliard.C
  LIBRARIES
    CommonFrac
)

ifem_add_application(
  NAME
    FractureDynamics
  SOURCES
    FractureArgs.C
    main_FractureDynamics.C
  HEADERS
    FractureArgs.h
  LIBRARIES
    CommonFrac
)

# Testing
enable_testing()

# Regression tests
ifem_add_hdf5_test(
  TARGET
    FractureDynamics
  TEST_FILES
    Sneddon-200_FD.hreg
)

ifem_add_hdf5_test(
  TARGET
    FractureDynamics
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Short10x20-p1-adap.hreg
)

ifem_add_vtf_test(
  TARGET
    FractureDynamics
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Short10x20-p1-adap.vreg
)

ifem_add_restart_test(
  TARGET
    FractureDynamics
  CONDITIONS
    LRSpline_FOUND
    Zoltan_FOUND
  TEST_FILES
    MPI/Short10x20-p1-adap-restart.reg
  RESTART_LEVEL
    38
  PARALLEL
    4
)

ifem_add_restart_test(
  TARGET
    FractureDynamics
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Short10x20-p1-adap-restart.reg
  RESTART_LEVEL
    38
)

ifem_add_restart_test(
  TARGET
    FractureDynamics
  TEST_FILES
    Miehe71_restart.reg
  RESTART_LEVEL
    11
)

ifem_add_restart_test(
  TARGET
    FractureDynamics
  CONDITIONS
    PoroElastic_DIRECTORY
  TEST_FILES
    Miehe71_restart_poro.reg
  RESTART_LEVEL
    11
)

ifem_add_regression_test(
  TARGET
    FractureDynamics
  TEST_FILES
    Miehe71-explcrack_FD.reg
    Miehe71_FD.reg
    Miehe71-lstatic_FD.reg
    Miehe71-nocrack_FD.reg
    Short10x20_FD.reg
    Sneddon-200_FD.reg
    ts1-25x25-5P-adap2-relax-aitken_FD.reg
    ts1-25x25-5P-adap2-relax_FD.reg
    ts1-50x40-Miehe_FD.reg
)

ifem_add_regression_test(
  TARGET
    FractureDynamics
  CONDITIONS
      LRSpline_FOUND
  TEST_FILES
    Short10x20-adap_FD.reg
    Short10x20-preref_FD.reg
)

ifem_add_regression_test(
  TARGET
    FractureDynamics
  CONDITIONS
    PoroElastic_DIRECTORY
  TEST_FILES
    # TODO: Re-enable these tests when poroelastic quasistatic is sensible
    #Miehe71-explcrack_FDporo.reg
    #Miehe71_FDporo.reg
    Miehe71-lstatic_FDporo.reg
    Miehe71-nocrack_FDporo.reg
    Sneddon-unit_FDporo.reg
)

# Cahn-Hilliard
ifem_add_hdf5_test(
  TARGET
    CahnHilliard
  TEST_FILES
    Rectangle-p2_CH.hreg
)

ifem_add_regression_test(
  TARGET
    CahnHilliard
  TEST_FILES
    Rectangle-p2_CH.reg
    Rectangle-p2_penalty_CH.reg
)

# Unit tests
ifem_add_test_app(
  NAME
    FractureDynamics
  SOURCES
    Test/TestFractureElasticity.C
    Test/TestSIMDynElasticity.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/Test
  LIBRARIES
    CommonFrac
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(
  TARGET
    OpenFrac
  DOX
    FractureDynamics
  EXTRA_DIRS
    ${Elasticity_DIRECTORY}
    ${PoroElastic_DIRECTORY}
)

# Installation
include(GNUInstallDirs)
install(TARGETS CahnHilliard FractureDynamics DESTINATION ${CMAKE_INSTALL_BINDIR})
