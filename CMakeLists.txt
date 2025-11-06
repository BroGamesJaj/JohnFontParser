cmake_minimum_required(VERSION 3.20)
project(JohnFont C)

set(CMAKE_C_STANDARD_REQUIRED ON)

file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/Source/*.c")

set(outputDir ${CMAKE_SOURCE_DIR}/Output)
if(NOT MSVC)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(outputDir ${CMAKE_SOURCE_DIR}/Output/Release)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(outputDir ${CMAKE_SOURCE_DIR}/Output/Debug)
    endif()
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${outputDir})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${outputDir})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${outputDir})

add_executable(JohnFont ${SOURCES})

target_include_directories(JohnFont PRIVATE ${CMAKE_SOURCE_DIR}/Source)

if(NOT MSVC)
    target_link_libraries(JohnFont PRIVATE m)  # link math library
endif()

