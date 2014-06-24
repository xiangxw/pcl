cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(slic_superpixel_segmentation)

find_package(PCL 1.7.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (slic_superpixel_segmentation slic_superpixel_segmentation.cpp)
target_link_libraries (slic_superpixel_segmentation ${PCL_LIBRARIES})
