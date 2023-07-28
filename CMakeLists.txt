cmake_minimum_required(VERSION 3.20)

## make sure to use --triplet=x64-windows when doing vcpckg install

set(CMAKE_BUILD_TYPE "Release")
# set the project name
project(Analysis)

set(OpenCV_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/opencv4_x64-windows/share/opencv/")
#set(unofficial-graphicsmagick_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/graphicsmagick_x86-windows/share/unofficial-graphicsmagick")


#set(LIBTIFF_INCL "${PROJECT_SOURCE_DIR}/libtiff/libtiff/")

#add_subdirectory(libtiff)

# test2
find_library(TIFF "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/bin/tiff.dll")


#find_library(BZIP2 "${PROJECT_SOURCE_DIR}/vcpkg/packages/bzip2_x86-windows/bin/bz2.dll")

# find_library(graphicsmagick "${PROJECT_SOURCE_DIR}/vcpkg/packages/graphicsmagick_x86-windows/bin/graphicsmagick.dll")
#set(TIFF_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/include/")
#set(TIFF_LIBRARY "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/lib/tiffxx.lib")


#set(BZIP2_LIBRARIES "${PROJECT_SOURCE_DIR}/vcpkg/packages/bzip2_x86-windows/lib/bz2.lib")
#set(BZIP2_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/bzip2_x86-windows/include/")
find_package(unofficial-graphicsmagick CONFIG REQUIRED)
find_package(OpenCV CONFIG REQUIRED)
find_package(TIFF REQUIRED) 

#
#add_executable(Analysis analysis.cpp)

set(EXECUTABLE_LIST analysis_endogenous_0x analysis2015-05-18 analysis_0x analysisU2OSHCR analysis_1x_10x tests mousebrain0207 mousebrain0822 example analysisU2OSHCR2 analysisU2OSHCR2023-07-11)

FOREACH(EXECUTABLE ${EXECUTABLE_LIST})
	add_executable(${EXECUTABLE} ${EXECUTABLE}.cpp)

#
	target_link_libraries(${EXECUTABLE} PRIVATE
					${OpenCV_LIBS}
					${TIFF_LIBS}
					${TIFF_LIBRARY}
					${unofficial-graphicsmagick_LIBS}
	)
	target_link_libraries(${EXECUTABLE} PUBLIC
					"${PROJECT_SOURCE_DIR}/LimFileStatic.lib"
					"${PROJECT_SOURCE_DIR}/Nd2ReadSdkStatic.lib"
					"${PROJECT_SOURCE_DIR}/LimFile.lib"
					"${PROJECT_SOURCE_DIR}/Nd2ReadSdk.lib"
					)

	target_include_directories(${EXECUTABLE} PUBLIC
						${PROJECT_BINARY_DIR}
						${OpenCV_INCLUDE_DIRS}
						${TIFF_INCLUDE_DIRS}
						${unofficial-graphicsmagick_INCLUDE_DIRS}
						) 

	set_property (TARGET ${EXECUTABLE}
		PROPERTY
		# Enable C++17 standard compliance
			CXX_STANDARD 17
	)
ENDFOREACH()
