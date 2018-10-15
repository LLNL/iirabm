
cmake_minimum_required( VERSION 2.8 )

project( BoostPythonSepsis )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -march=native")

# Mac OSX-specific fixes
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
  set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif()

# Find necessary packages
find_package( PythonLibs 3.6 REQUIRED )
include_directories( ${PYTHON_INCLUDE_DIRS} )

find_package( Boost COMPONENTS python3 numpy3 REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

file(GLOB_RECURSE SepsisSRC "src/core/*.cpp")
file(GLOB_RECURSE SepsisHEADERS "src/core/*.h")

set (Sepsis_INCLUDE_DIRS "")
foreach (_headerFile ${SepsisHEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND Sepsis_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES Sepsis_INCLUDE_DIRS)

# Build our library
add_library( sepsis SHARED src/simulationControl.cpp ${SepsisSRC} )
target_include_directories(sepsis PRIVATE ${Sepsis_INCLUDE_DIRS})

target_link_libraries( sepsis ${Boost_LIBRARIES} )
set_target_properties( sepsis PROPERTIES PREFIX "" )
