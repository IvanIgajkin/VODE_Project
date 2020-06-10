cmake_minimum_required(VERSION 3.10)
project(ODE)

enable_language(Fortran)

set(CMAKE_MAKE_PROGRAM 'make')

add_executable(ODE
        source.cpp
        vode.f
        ode.f90
        linpack/dgbsl.f
        linpack/dgesl.f
        linpack/dgefa.f
        linpack/dgbfa.f)

target_link_libraries(ODE blas)
