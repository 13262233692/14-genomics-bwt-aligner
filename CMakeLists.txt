cmake_minimum_required(VERSION 3.15)
project(genomics_bwt_aligner VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    add_compile_options(/O2 /W3 /MD)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        add_compile_options(/arch:AVX2)
    endif()
else()
    add_compile_options(-O3 -Wall -Wextra -pthread)
    add_link_options(-pthread)
    add_compile_options(-mavx2 -mfma)
endif()

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CORE_SOURCES
    ${SRC_DIR}/fasta_reader.cpp
    ${SRC_DIR}/fastq_reader.cpp
    ${SRC_DIR}/suffix_array.cpp
    ${SRC_DIR}/bwt.cpp
    ${SRC_DIR}/fm_index.cpp
    ${SRC_DIR}/backward_search.cpp
    ${SRC_DIR}/sam_writer.cpp
    ${SRC_DIR}/smith_waterman.cpp
    ${SRC_DIR}/aligner.cpp
)

set(CORE_HEADERS
    ${INCLUDE_DIR}/fasta_reader.h
    ${INCLUDE_DIR}/fastq_reader.h
    ${INCLUDE_DIR}/suffix_array.h
    ${INCLUDE_DIR}/bwt.h
    ${INCLUDE_DIR}/fm_index.h
    ${INCLUDE_DIR}/backward_search.h
    ${INCLUDE_DIR}/sam_writer.h
    ${INCLUDE_DIR}/smith_waterman.h
    ${INCLUDE_DIR}/aligner.h
    ${INCLUDE_DIR}/sequence.h
)

add_library(bwt_aligner_core STATIC ${CORE_SOURCES} ${CORE_HEADERS})
target_include_directories(bwt_aligner_core PUBLIC ${INCLUDE_DIR})

option(BUILD_PYTHON_BINDINGS "Build Python bindings with pybind11" ON)

if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter Development)
    if(Python3_FOUND)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
            OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE PYBIND11_RESULT
        )
        if(NOT PYBIND11_RESULT EQUAL 0)
            message(WARNING "pybind11 not found via Python, trying find_package")
            find_package(pybind11 CONFIG)
        else()
            list(APPEND CMAKE_PREFIX_PATH ${PYBIND11_CMAKE_DIR})
            find_package(pybind11 CONFIG REQUIRED)
        endif()

        if(pybind11_FOUND)
            pybind11_add_module(bwt_aligner ${SRC_DIR}/pybind_module.cpp)
            target_link_libraries(bwt_aligner PRIVATE bwt_aligner_core)
            target_include_directories(bwt_aligner PRIVATE ${INCLUDE_DIR})
        else()
            message(WARNING "pybind11 not found, skipping Python bindings")
        endif()
    else()
        message(WARNING "Python not found, skipping Python bindings")
    endif()
endif()

option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

add_executable(test_sa test_sa.cpp)
target_link_libraries(test_sa bwt_aligner_core)
target_include_directories(test_sa PRIVATE ${INCLUDE_DIR})

add_executable(test_sw test_sw.cpp)
target_link_libraries(test_sw bwt_aligner_core)
target_include_directories(test_sw PRIVATE ${INCLUDE_DIR})
