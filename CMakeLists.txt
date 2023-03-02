cmake_minimum_required(VERSION 3.20)

set(CMAKE_BUILD_TYPE "Release")
# set the project name
project(Analysis)

set(OpenCV_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/opencv4_x64-windows/share/opencv/")
#set(LIBTIFF_INCL "${PROJECT_SOURCE_DIR}/libtiff/libtiff/")

#add_subdirectory(libtiff)

## test
find_library(TIFF "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/bin/tiff.dll")
#set(TIFF_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/include/")
#set(TIFF_LIBRARY "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/lib/tiffxx.lib")


find_package(OpenCV CONFIG REQUIRED)
find_package(TIFF REQUIRED) 

#
#add_executable(Analysis analysis.cpp)

set(EXECUTABLE_LIST tests example)

FOREACH(EXECUTABLE ${EXECUTABLE_LIST})
	add_executable(${EXECUTABLE} ${EXECUTABLE}.cpp)

#
	target_link_libraries(${EXECUTABLE} PRIVATE
					${OpenCV_LIBS}
					${TIFF_LIBS}
					${TIFF_LIBRARY}
					)

	target_include_directories(${EXECUTABLE} PUBLIC
						${PROJECT_BINARY_DIR}
						${OpenCV_INCLUDE_DIRS}
						${TIFF_INCLUDE_DIRS}
						) 

	set_property (TARGET ${EXECUTABLE}
		PROPERTY
		# Enable C++17 standard compliance
			CXX_STANDARD 17
	)
ENDFOREACH()
