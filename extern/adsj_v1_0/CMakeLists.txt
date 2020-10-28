cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(adsj VERSION 1.0.0 LANGUAGES Fortran C CXX)

include(FortranCInterface)
FortranCInterface_VERIFY()

if (NOT CMAKE_BUILD_TYPE)
    # other options: Debug, Release, RelWithDebInfo, MinSizeRel
    #set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-fno-range-check" has_fno_range_check)
if(has_fno_range_check)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()
check_fortran_compiler_flag("-fallow-argument-mismatch" has_fallow_argument_mismatch)
if(has_fallow_argument_mismatch)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif()

add_library(adsj STATIC
    check_parms.f90
    driver.f90
    f1bnl.f90
    f1conv.f90
    f1v3.f90
    f2.f90
    f21r.f
    f2gam.f90
    gamma.f90
    get_transformation2.f90
    hypdrvf1.f90
    hypgeof1.f90
    hypstartf1.f90
    isneg.f90
    writef1.f90
    webroutines/chyp.f
    webroutines/d1mach.for
    webroutines/f1dxy.f90
    webroutines/hyp.f
    webroutines/hypgeo.f
    webroutines/rkf45.for
    adsj_c_wrapper.f90
)

add_executable(driver driver.f90)

target_link_libraries(driver PUBLIC
    adsj
)
