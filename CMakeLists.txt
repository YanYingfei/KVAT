cmake_minimum_required (VERSION 3.18)

project (Anonymous-Tokens-with-keyed-verification LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectory(poly)
add_subdirectory(test)
add_subdirectory(protocol)