
cmake_minimum_required(VERSION "3.14")

project("rad_hydro")

#set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-10" CACHE STRING "C++ compiler" FORCE)

set(SOURCE_FILES "src/main.cpp")

#add_compile_options(-Wextra -Wall -Wpedantic)

#Required for Lua to run
find_package(Lua51 REQUIRED)
include_directories(/usr/local/include ${LUA_INCLUDE_DIR})

add_executable("${PROJECT_NAME}" "${SOURCE_FILES}")
target_link_libraries("${PROJECT_NAME}" ${LUA_LIBRARIES})
