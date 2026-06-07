cmake_minimum_required(VERSION 3.16)
project(EigenHDF5 LANGUAGES C CXX CUDA)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -allow-unsupported-compiler")

# Tell nvcc which host compiler to use — must match what compiled st1.cpp
set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

add_compile_definitions(EIGEN_DONT_VECTORIZE)
add_compile_definitions(EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT)

# --- MKL Configuration ---
# set(MKLROOT "/opt/intel/oneapi/mkl/2026.0")
# set(INTEL_LIB_PATH "/opt/intel/oneapi/compiler/latest/lib")
set(MKLROOT "/apps/compilers/intel/oneapi/2023.2.0/mkl/2023.2.0")
set(INTEL_LIB_PATH "/apps/compilers/intel/oneapi/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin")
set(CUDAToolkit_ROOT "/usr/local/cuda/bin")

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(CUDAToolkit REQUIRED)

add_executable(st st.cpp Master/Operators/GPU.cu Master/Operators/GPU.cpp)
add_executable(mt mt.cpp Master/Operators/GPU.cu Master/Operators/GPU.cpp)

# Common include directories
set(COMMON_INCLUDES
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/Master
    ${CMAKE_SOURCE_DIR}/Master/Eigen
    ${CMAKE_SOURCE_DIR}/Master/HighFive/include
    ${HDF5_INCLUDE_DIRS}
    ${MKLROOT}/include
)

target_include_directories(st PRIVATE ${COMMON_INCLUDES})
target_include_directories(mt PRIVATE ${COMMON_INCLUDES})

# Link st with MKL and OpenMP
target_link_libraries(st PRIVATE
    Threads::Threads
    HDF5::HDF5
    OpenMP::OpenMP_CXX
    -L${MKLROOT}/lib/intel64
    -L${INTEL_LIB_PATH}
    -Wl,--start-group
      ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
      ${MKLROOT}/lib/intel64/libmkl_core.a
      ${MKLROOT}/lib/intel64/libmkl_intel_thread.a
    -Wl,--end-group
    iomp5
    pthread
    m
    dl
    CUDA::cudart
    CUDA::cublas
    CUDA::cusolver
)

target_link_libraries(mt PRIVATE
    Threads::Threads
    HDF5::HDF5
    OpenMP::OpenMP_CXX
    -L${MKLROOT}/lib/intel64
    -L${INTEL_LIB_PATH}
    -Wl,--start-group
      ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
      ${MKLROOT}/lib/intel64/libmkl_core.a
      ${MKLROOT}/lib/intel64/libmkl_intel_thread.a
    -Wl,--end-group
    iomp5
    pthread
    m
    dl
    CUDA::cudart
    CUDA::cublas
    CUDA::cusolver
)

set_target_properties(st PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES 89)

set_target_properties(mt PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES 89)
