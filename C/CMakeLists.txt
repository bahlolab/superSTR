cmake_minimum_required(VERSION 3.10)    # CMake version check
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Cmakes)

add_definitions(-g3 -g)
#add_definitions(-fsanitize=address -g3 -fno-omit-frame-pointer)
#add_compile_options(-ggdb -fsanitize=address)
#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

project(superSTR)                 # Create project "simple_example"
set(CMAKE_C_STANDARD 99)                # Enable C99 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.c factorizeforGDR.c finalstageforGDR.c findReps.c kmp.c mainRepeats.c mainsearchforGDR.c mreps.c printOutput.c searchforHeadGDR.c argtable3.c)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(superstr ${SOURCE_FILES})

find_package(ZLIB)

if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries(superstr ${ZLIB_LIBRARIES})
endif()

find_package(CURL)

if (CURL_FOUND)
    include_directories(${CURL_INCLUDE_DIRS})
    target_link_libraries(superstr ${CURL_LIBRARIES})
endif()

find_package(htslib)

if (HTSlib_FOUND)
    include_directories(${HTSlib_INCLUDE_DIRS})
    target_link_libraries(superstr ${HTSlib_LIBRARIES})
endif()

target_link_libraries(superstr m)
