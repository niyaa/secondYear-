CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# TODO - Put the name of your project here.
SET(ProjectName Nikesh)
PROJECT(${ProjectName})
# Setup where CMake will look for modules.  You probably
# won't need to modify this.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})
SET(INCLUDE_DIRS "/usr/include/eigen3/Eigen")
ADD_SUBDIRECTORY(source1d)

