cmake_minimum_required(VERSION 2.8)

project(fkrig)

# add_definitions("-std=c++11 -Wall -Wunused-private-field -g -O0")
add_definitions("-std=c++11 -Wall -g -Wno-sign-compare -O0")

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

include_directories(/usr/local/include/mba)

include_directories(./include)

add_library( fkrig ./src/util_fkrig.cpp ./src/covariance.cpp ./src/curve_base.cpp ./src/curve.cpp ./src/surface_base.cpp ./src/surface.cpp ./src/ego_base.cpp ./src/ego_curve.cpp ./src/ego_surface.cpp ./src/curve_link.cpp ./src/surface_link.cpp)

# target_link_libraries(fkrig ${CERES_LIBRARIES})

# add_executable(fkrig main.cpp)

# install(TARGETS fkrig RUNTIME DESTINATION bin)