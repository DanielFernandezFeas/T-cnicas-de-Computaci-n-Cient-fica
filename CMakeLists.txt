cmake_minimum_required(VERSION 3.20)
project(heat2d C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Carpeta de headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Código común (malla, IO, timer...)
add_library(common
        src/grid.c
        src/io.c
        src/timer.c
        "Proyecto Secuencial/include/config.h"
        "Proyecto Secuencial/include/config.h"
        "Proyecto Secuencial/include/grid.h"
        "Proyecto Secuencial/src/heat_seq_opt.c"
)

# ---------- SEQ ----------
add_executable(main_seq
        apps/main_seq.c
        src/heat_seq.c
)
target_link_libraries(main_seq PRIVATE common)

# ---------- OpenMP ----------
find_package(OpenMP)
if(OpenMP_C_FOUND)
    add_executable(main_omp
            apps/main_omp.c
            src/heat_omp.c
    )
    target_link_libraries(main_omp PRIVATE common OpenMP::OpenMP_C)
else()
    message(WARNING "OpenMP no encontrado: se omitirá main_omp")
endif()

# ---------- MPI ----------
find_package(MPI)
if(MPI_C_FOUND)
    add_executable(main_mpi
            apps/main_mpi.c
            src/heat_mpi.c
    )
    target_link_libraries(main_mpi PRIVATE common MPI::MPI_C)
else()
    message(WARNING "MPI no encontrado: se omitirá main_mpi")
endif()
