#include <iostream>
#include <cmath>
#include <tuple>
#include <optional>
#include <stack>
#include <concurrent_vector.h>
#include <time.h>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/core/mat.hpp>
#include <tiff.h>
#include <tiffio.h>
#include <fstream>
// WINDOWS SPECIFIC CONCURRENCY RUNTIME
#include <ppl.h>
#include <ppltasks.h>


#define M_PI 3.14159265358979323846

/*  
   Painful process learning how to use CMAKE, but vcpkg seems to make it very easy to install new ones, just use ./vcpkg/vcpkg install. 

   In order to do this again:
   - install packages with vcpkg
   - set up CMakeLists.txt file w/ calls to find_package and target_link_libraries

   DO NOT INSTALL/BUILD INSIDE EMACS. USE FULL SHELL. OR SUFFER PAIN. 

   To build:
   cd build
   cmake .. -DCMAKE_TOOLCHAIN_FILE=/home/andreyshur/mike_cant_access_his_externaldrive/ImgAnalysis/vcpkg/scripts/buildsystems/vcpkg.cmake 
   make

   To set up with flycheck, use -DCMAKE_EXPORT_COMPILE_COMMANDS=ON and look in compile_commands.json, put whatever is -isystem into .dir-locals.el flycheck-gcc-include-path

   GUI TOOL INSTALLATION IS HELL 
   Installing gtk with vcpkg is hell. When things fail it's because
   libraries are missing, like libxi-dev and other, and must be
   install with apt.

   For Qt, a list of system requirements to be installed is here:
   https://doc.qt.io/qt-5/linux-requirements.html

   Might need this for the many xcb dependencies
   sudo apt-get install '^libxcb.*-dev' libx11-xcb-dev libglu1-mesa-dev libxrender-dev libxi-dev libxkbcommon-dev libxkbcommon-x11-dev

 */

#define GAUSSIAN_CUTOFF 4
#define zincrement 2048*2048
#define yincrement 2048

struct BGR {
	uint16_t blue;
	uint16_t green;
	uint16_t red;
};

struct BGR_float {
	float blue;
	float green;
	float red;
};

struct float3 {
	float x;
	float y;
	float z; 
};

class Blob {
public:
	std::vector<std::tuple<int, int, int>> points;
	std::vector<std::tuple<int, int, int>> boundary;
	std::tuple<int, int, int> local_max;

	int size() {
		return points.size();
	}
	
	//virtual void to_csv(char buf[]) = 0;
};
class Dot;
class Nucleus : public Blob {
	const int size_upper_limit = 1000000;
	const int size_lower_limit = 100000;
public:
	int id;
	std::vector<Dot*> close_dots594;
	std::vector<Dot*> close_dots640;

	bool validSize() {
		return points.size() > size_lower_limit && points.size() < size_upper_limit;
	}
};
class Dot : public Blob {
	const int size_upper_limit = 5000;
	const int size_lower_limit = 1;
public:
	int id; 
	
	bool validSize() {
		return points.size() > size_lower_limit && points.size() < size_upper_limit;
	}
};

void gaussianElimination(std::array<std::array<float, 3>,3> mat) {
	// this typechecks?  
	return;
}

void findMaxima(float* voxels, std::vector<std::tuple<int, int, int>> &maxima) {
	for (int z = 2; z < 199; z++) {
		for (int y = 2; y < 2046; y++) {
			for (int x = 2; x < 2046; x++) {
				float current = voxels[2048 * y + 2048 * 2048 * z + x];
				// 6 checks
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + yincrement + zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + yincrement - zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - yincrement + zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - yincrement - zincrement])
					continue;

				// go out to distance of 2 to avoid numerical error

				// x - 2 face
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 2*yincrement - 2*zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 2 * yincrement ])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 2 * yincrement + 2* zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 1 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 1 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 1 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 1 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 - 1 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2  - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2  - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 ])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2  + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 1 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 1 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 1 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 1 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 +  1 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 + 2 * yincrement + 2 * zincrement])
					continue;

				// x + 2 face
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 2 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 1 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 1 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 1 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 1 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 1 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 1 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 1 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 1 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 1 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 1 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 + 2 * yincrement + 2 * zincrement])
					continue;


				// y - 2 face
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - 2 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x  - 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x  - 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x  - 2 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 2 * yincrement + 2 * zincrement])
					continue;


				// y + 2 face
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 2 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 2 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 2 * yincrement - 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 2 * yincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 2 * yincrement + 1 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 2 * yincrement + 2 * zincrement])
					continue;

				// z - 2 face
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - 1 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1  - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 1 * yincrement - 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 * yincrement - 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 1 * yincrement - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 1 * yincrement - 2 * zincrement])
					continue;


				// z + 2 face
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 - 1 * yincrement + 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 + 1 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x - 1 * yincrement + 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 * yincrement + 2 * zincrement])
					continue;

				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 - 1 * yincrement + 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 2 * zincrement])
					continue;
				if (current <= voxels[2048 * y + 2048 * 2048 * z + x + 1 + 1 * yincrement + 2 * zincrement])
					continue;

				maxima.push_back(std::make_tuple(x, y, z));
			}
		}
	}
}

void gradientField3d(float* voxels, float3* gradientField) {
	float sixth_order_centered[4] = {0, 3.0 / 4, -3.0 / 20, 1.0 / 60};
	float sixth_order_forward[7] = {-49.0 / 20, 6.0, -15.0 / 2, 20.0 / 3, -15.0 / 4, 6.0 / 5, -1.0 / 6};
	concurrency::parallel_for(0, 201, [&sixth_order_centered, &sixth_order_forward, gradientField, voxels](int z) {
		for (int y = 0; y < 2048; y++) {
			for (int x = 0; x < 2048; x++) {
				float dx = 0, dy = 0, dz = 0;

				if (x < 3) {
					for (int i = 0; i < 7; i++) {
						dx += voxels[(x + i) + 2048 * y + 2048 * 2048 * z] * sixth_order_forward[i];
					}
				}
				else if (x > 2048 - 4) {
					for (int i = 0; i < 7; i++) {
						dx += -voxels[(x - i) + 2048 * y + 2048 * 2048 * z] * sixth_order_forward[i];
					}
				}
				else {
					for (int i = 1; i < 4; i++) {
						dx += voxels[(x + i) + 2048 * y + 2048 * 2048 * z] * sixth_order_centered[i];
						dx += -voxels[(x - i) + 2048 * y + 2048 * 2048 * z] * sixth_order_centered[i];
					}
				}
				gradientField[x + 2048 * y + 2048 * 2048 * z].x = dx;

				if (y < 3) {
					for (int j = 0; j < 7; j++) {
						dy += voxels[x + 2048 * (y + j) + 2048 * 2048 * z] * sixth_order_forward[j];
					}
				}
				else if (y > 2048 - 4) {
					for (int j = 0; j < 7; j++) {
						dy += -voxels[x + 2048 * (y - j) + 2048 * 2048 * z] * sixth_order_forward[j];
					}
				}
				else {
					for (int j = 1; j < 4; j++) {
						dy += voxels[x + 2048 * (y + j) + 2048 * 2048 * z] * sixth_order_centered[j];
						dy += -voxels[x + 2048 * (y - j) + 2048 * 2048 * z] * sixth_order_centered[j];
					}
				}
				gradientField[x + 2048 * y + 2048 * 2048 * z].y = dy;

				if (z < 3) {
					for (int k = 0; k < 7; k++) {
						dz += voxels[x + 2048 * y + 2048 * 2048 * (z + k)] * sixth_order_forward[k];
					}
				}
				else if (z > 201 - 4) {
					for (int k = 0; k < 7; k++) {
						dz += -voxels[x + 2048 * y + 2048 * 2048 * (z - k)] * sixth_order_forward[k];
					}
				}
				else {
					for (int k = 1; k < 4; k++) {
						dz += voxels[x + 2048 * y + 2048 * 2048 * (z + k)] * sixth_order_centered[k];
						dz += -voxels[x + 2048 * y + 2048 * 2048 * (z - k)] * sixth_order_centered[k];

					}
				}
				gradientField[x + 2048 * y + 2048 * 2048 * z].z = dz;
			}
		}
		});
}

void eigenvalues(std::array<std::array<float, 3>, 3> hessian, float& eig1, float& eig2, float&eig3) {
	/*
*
*
	Extra Fast algorithm for computing eigenvalues of 3x3 symmetric matrices
	https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
*/
	float q = (hessian[0][0] + hessian[1][1] + hessian[2][2]) / 3;
	float p2 = (hessian[0][0] - q) * (hessian[0][0] - q) +
		(hessian[1][1] - q) * (hessian[1][1] - q) +
		(hessian[2][2] - q) * (hessian[2][2] - q) +
		2 * (hessian[1][0] * hessian[1][0] + hessian[2][0] * hessian[2][0] + hessian[2][1] * hessian[2][1]);
	float p = sqrt(p2 / 6);

	// B = 1/p (A - qI)
	hessian[0][0] -= q;
	hessian[1][1] -= q;
	hessian[2][2] -= q;

	hessian[0][0] = hessian[0][0] / p;
	hessian[1][1] = hessian[1][1] / p;
	hessian[2][2] = hessian[2][2] / p;
	hessian[1][0] = hessian[1][0] / p;
	hessian[2][0] = hessian[2][0] / p;
	hessian[2][1] = hessian[2][1] / p;


	float det = hessian[0][0] * (hessian[1][1] * hessian[2][2] - hessian[2][1] * hessian[2][1])
		- hessian[1][0] * (hessian[1][0] * hessian[2][2] - hessian[2][1] * hessian[2][0])
		+ hessian[2][0] * (hessian[1][0] * hessian[2][1] - hessian[1][1] * hessian[2][0]);


	float r2 = det / 2;
	float phi;
	if (r2 <= -1) {
		phi = 3.1415926535 / 3;
	}
	else if (r2 >= 1.0) {
		phi = 0;
	}
	else {
		phi = acos(r2) / 3;
	}
	eig1 = q + 2 * p * cos(phi);
	eig2 = q + 2 * p * cos(phi + (2 * 3.1415926535 / 3));
	eig3 = 3 * q - eig1 - eig2;
}

// this is not the correct way to compute the eigenvector, but it'll have to do for now.
// why not just swap this out for more professional algos like LAPACK?
void eigenvector(std::array<std::array<float, 3>, 3> hessian, float eig, float3& ev) {
	hessian[0][0] = hessian[0][0] - eig;
	hessian[1][1] = hessian[1][1] - eig;
	hessian[2][2] = hessian[2][2] - eig;

	std::array<float, 3> tmp{};

	for (int i = 0; i < 3; i++) {
		// change the pivot if the current point is zero
		if(hessian[i][i] == 0) {
			for (int j = i+1; j < 3; j++) {
				// swap pivot row into ith row
				if (hessian[j][i] != 0) {
					for (int k = 0; k < 3; k++) {
						tmp[k] = hessian[i][k];
						hessian[i][k] = hessian[j][k]; // normalize to 1
						hessian[j][k] = tmp[k];
					}
					break;
				}
			}  // if we get to the end of this, it's just a column of zeros.
		}
		if (hessian[i][i] != 0) {
			for (int j = 0; j < 3; j++) {
				float ratio = hessian[j][i] / hessian[i][i];
				if (j == i) continue;
				for (int k = 0; k < 3; k++) {
					hessian[j][k] -= ratio * hessian[i][k];
				} // zeroing out the jth column
			}
		}


	} // now we should be in some kind of rref


	for (int i = 0; i < 3; i++) {
		if (hessian[i][i] == 0) {
			tmp[i] = 1;
		}
		else {
			for (int j = 0; j < 3; j++) {
				if(j!=i) tmp[i] -= hessian[i][j] / hessian[i][i];
			}
		}
	}

	float ev_mag = sqrt(tmp[0]*tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);

	ev.x = tmp[0] / ev_mag;
	ev.y = tmp[1] / ev_mag;
	ev.z = tmp[2] / ev_mag;
}


// does not compute laplacian on the boundary.
// can be implemented later with a one-sided stencil. 
// but I could not find the terms using a quick google search. 
void laplacianFilter3D(float* voxels, float* laplacian) {
	// z x y laplacian
	float sixth_order_centered[4] = {-49.0 / 18, 3.0 / 2, -3.0 / 20, 1.0 / 90};
	float sixth_order_forward[8] = {469.0 / 90, -223.0 / 10, 879.0 / 20, -949.0 / 18, 41.0, -201.0 / 10, 1019.0 / 180, -7.0 / 10};

	for (int z = 0; z < 201; z++) {
		for (int y = 0; y < 2048; y++) {
			for (int x = 0; x < 2048; x++) {
				float lp = 0;
				// x
				if (x < 3) {
					for (int i = 0; i < 8; i++) {
						lp += voxels[(x + i) + 2048 * y + 2048 * 2048 * z] * sixth_order_forward[i];
					}
				}
				else if (x > 2048 - 4) {
					for (int i = 0; i < 8; i++) {
						lp += voxels[(x - i) + 2048 * y + 2048 * 2048 * z] * sixth_order_forward[i];
					}
				}
				else {
					lp += voxels[x + 2048 * y + 2048 * 2048 * z] * sixth_order_centered[0];
					for (int i = 1; i < 4; i++) {
						lp += voxels[(x+i) + 2048 * y + 2048 * 2048 * z] * sixth_order_centered[i];
						lp += voxels[(x-i) + 2048 * y + 2048 * 2048 * z] * sixth_order_centered[i];
					}
				}

				// y
				if (y < 3) {
					for (int j = 0; j < 8; j++) {
						lp += voxels[x + 2048 * (y+j) + 2048 * 2048 * z] * sixth_order_forward[j];
					}
				} 
				else if (y > 2048 - 4) {
					for (int j = 0; j < 8; j++) {
						lp += voxels[x + 2048 * (y - j) + 2048 * 2048 * z] * sixth_order_forward[j];
					}
				}
				else {
					lp += voxels[x + 2048 * y + 2048 * 2048 * z] * sixth_order_centered[0];
					for (int j = 1; j < 4; j++) {
						lp += voxels[x + 2048 * (y + j) + 2048 * 2048 * z] * sixth_order_centered[j];
						lp += voxels[x + 2048 * (y - j) + 2048 * 2048 * z] * sixth_order_centered[j];

					}
				}

				// z 
				if (z < 3) {
					for (int k = 0; k < 8; k++) {
						lp += voxels[x + 2048 * y + 2048 * 2048 * (z+k)] * sixth_order_forward[k];

					}
				}
				else if (z > 201 - 4) {
					for (int k = 0; k < 8; k++) {
						lp += voxels[x + 2048 * y + 2048 * 2048 * (z - k)] * sixth_order_forward[k];

					}
				}
				else {
					lp += voxels[x + 2048 * y + 2048 * 2048 * z] * sixth_order_centered[0];
					for (int k = 1; k < 4; k++) {
						lp += voxels[x + 2048 * y + 2048 * 2048 * (z+k)] * sixth_order_centered[k];
						lp += voxels[x + 2048 * y + 2048 * 2048 * (z-k)] * sixth_order_centered[k];
					}
				}
				laplacian[x + 2048 * y + 2048 * 2048 * z] = lp;
			}
		}
	}
}

void hessianAt(float3* gradientField, int x2, int y2, int z2, std::array<std::array<float, 3>, 3>& hessian) {
	float sixth_order_centered[4] = { 0, 3.0 / 4, -3.0 / 20, 1.0 / 60 };
	float sixth_order_forward[7] = { -49.0 / 20, 6.0, -15.0 / 2, 20.0 / 3, -15.0 / 4, 6.0 / 5, -1.0 / 6 };

	if (x2 < 3) {
		for (int j = 0; j < 7; j++) {
			float3 gradf = gradientField[(x2 + j) + 2048 * y2 + 2048 * 2048 * z2];
			hessian[0][0] += gradf.x * sixth_order_forward[j];// d2fdx2
			hessian[1][0] += gradf.y * sixth_order_forward[j]; // d2fdxdy
			hessian[2][0] += gradf.z * sixth_order_forward[j]; // d2fdxdz
		}
	}
	else if (x2 > 2048 - 4) {
		for (int j = 0; j < 7; j++) {
			float3 gradf = gradientField[(x2 - j) + 2048 * y2 + 2048 * 2048 * z2];
			hessian[0][0] -= gradf.x * sixth_order_forward[j];// d2fdx2
			hessian[1][0] -= gradf.y * sixth_order_forward[j]; // d2fdxdy
			hessian[2][0] -= gradf.z * sixth_order_forward[j]; // d2fdxdz
		}
	}
	else {
		for (int j = 1; j < 4; j++) {
			float3 gradf = gradientField[(x2 - j) + 2048 * y2 + 2048 * 2048 * z2];
			hessian[0][0] -= gradf.x * sixth_order_centered[j];// d2fdx2
			hessian[1][0] -= gradf.y * sixth_order_centered[j]; // d2fdxdy
			hessian[2][0] -= gradf.z * sixth_order_centered[j]; // d2fdxdz


			gradf = gradientField[(x2 + j) + 2048 * y2 + 2048 * 2048 * z2];
			hessian[0][0] += gradf.x * sixth_order_centered[j];// d2fdx2
			hessian[1][0] += gradf.y * sixth_order_centered[j]; // d2fdxdy
			hessian[2][0] += gradf.z * sixth_order_centered[j]; // d2fdxdz
		}
	}


	if (y2 < 3) {
		for (int k = 0; k < 7; k++) {
			float3 gradf = gradientField[x2 + 2048 * (y2 + k) + 2048 * 2048 * z2];
			hessian[1][1] += gradf.y * sixth_order_forward[k];// d2fdy2
			hessian[2][1] += gradf.z * sixth_order_forward[k];// d2fdydz
		}
	}
	else if (y2 > 2048 - 4) {
		for (int k = 0; k < 7; k++) {
			float3 gradf = gradientField[x2 + 2048 * (y2 - k) + 2048 * 2048 * z2];
			hessian[1][1] -= gradf.y * sixth_order_forward[k];// d2fdy2
			hessian[2][1] -= gradf.z * sixth_order_forward[k];// d2fdy2
		}
	}
	else {
		for (int k = 1; k < 4; k++) {
			float3 gradf = gradientField[x2 + 2048 * (y2 + k) + 2048 * 2048 * z2];
			hessian[1][1] += gradf.y * sixth_order_centered[k];// d2fdy2
			hessian[2][1] += gradf.z * sixth_order_centered[k];// d2fdy2


			gradf = gradientField[x2 + 2048 * (y2 - k) + 2048 * 2048 * z2];
			hessian[1][1] -= gradf.y * sixth_order_centered[k];// d2fdy2
			hessian[2][1] -= gradf.z * sixth_order_centered[k];// d2fdy2
		}
	}


	if (z2 < 3) {

		for (int l = 0; l < 7; l++) {
			float3 gradf = gradientField[x2 + 2048 * y2 + 2048 * 2048 * (z2 + l)];
			hessian[2][2] += gradf.z * sixth_order_forward[l];// d2fdz2

		}
	}
	else if (z2 > 201 - 4) {
		for (int l = 0; l < 7; l++) {
			float3 gradf = gradientField[x2 + 2048 * y2 + 2048 * 2048 * (z2 - l)];
			hessian[2][2] -= gradf.z * sixth_order_forward[l];// d2fdz2
		}
	}
	else {

		for (int l = 1; l < 4; l++) {
			float3 gradf = gradientField[x2 + 2048 * y2 + 2048 * 2048 * (z2 + l)];
			hessian[2][2] += gradf.z * sixth_order_centered[l];// d2fdz2

			gradf = gradientField[x2 + 2048 * y2 + 2048 * 2048 * (z2 - l)];
			hessian[2][2] -= gradf.z * sixth_order_centered[l];// d2fdz2

		}
	}
	hessian[0][1] = hessian[1][0];
	hessian[0][2] = hessian[2][0];
	hessian[1][2] = hessian[2][1];
}

// second order accurate
float third_directional_deriv(float3* gradientField, int x3, int y3, int z3, float3& ev) {
	// compute the hessian at x2
	float3 d3f{};
	std::array<std::array<float, 3>, 3> hessian_tmp;

	d3f.x = 0;
	if (x3 < 1) {
		hessianAt(gradientField, x3, y3, z3, hessian_tmp);
		d3f.x -= 3.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3 + 1, y3, z3, hessian_tmp);
		d3f.x += 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3 + 2, y3, z3, hessian_tmp);
		d3f.x -= 1.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	else if (x3 > 2048 - 2) {
		hessianAt(gradientField, x3, y3, z3, hessian_tmp);
		d3f.x += 3.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3 - 1, y3, z3, hessian_tmp);
		d3f.x -= 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3 - 2, y3, z3, hessian_tmp);
		d3f.x += 1.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	else {
		hessianAt(gradientField, x3 + 1, y3, z3, hessian_tmp);
		d3f.x += (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3 - 1, y3, z3, hessian_tmp);
		d3f.x -= (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}

	d3f.y = 0;
	if (y3 < 1) {
		hessianAt(gradientField, x3, y3, z3, hessian_tmp);
		d3f.y -= 3.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3 + 1, z3, hessian_tmp);
		d3f.y += 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3 + 2, z3, hessian_tmp);
		d3f.y -= 1.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	else if (y3 > 2048 - 2) {
		hessianAt(gradientField, x3, y3, z3, hessian_tmp);
		d3f.y += 3.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3 - 1, z3, hessian_tmp);
		d3f.y -= 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3 - 2, z3, hessian_tmp);
		d3f.y += 1.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	else {
		hessianAt(gradientField, x3, y3 + 1, z3, hessian_tmp);
		d3f.y += (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3 - 1, z3, hessian_tmp);
		d3f.y -= (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	d3f.z = 0;
	if (z3 < 1) {
		hessianAt(gradientField, x3, y3, z3, hessian_tmp);
		d3f.z -= 3.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3, z3 + 1, hessian_tmp);
		d3f.z += 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3, z3 + 2, hessian_tmp);
		d3f.z -= 1.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	else if (z3 > 201 - 2) {
		hessianAt(gradientField, x3, y3, z3, hessian_tmp);
		d3f.z += 3.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3, z3 - 1, hessian_tmp);
		d3f.z -= 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3, z3 - 2, hessian_tmp);
		d3f.z += 1.0 / 2 * (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	else {
		hessianAt(gradientField, x3, y3, z3 + 1, hessian_tmp);
		d3f.z += (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
		hessianAt(gradientField, x3, y3, z3 - 1, hessian_tmp);
		d3f.z -= (ev.x * (hessian_tmp[0][0] * ev.x + hessian_tmp[1][0] * ev.y + hessian_tmp[2][0] * ev.z) +
			ev.y * (hessian_tmp[1][0] * ev.x + hessian_tmp[1][1] * ev.y + hessian_tmp[2][1] * ev.z) +
			ev.z * (hessian_tmp[2][0] * ev.x + hessian_tmp[2][1] * ev.y + hessian_tmp[2][2] * ev.z));
	}
	return d3f.x * ev.x + d3f.y * ev.y + d3f.z * ev.z;
}

// second order accurate
float fourth_directional_deriv(float3* gradientField, int x, int y, int z, float3& ev) {
	float3 df4{};

	df4.x = 0;
	if (x < 1) {
		df4.x -= 3.0 / 2 * third_directional_deriv(gradientField, x, y, z, ev);
		df4.x += 2 * third_directional_deriv(gradientField, x + 1, y, z, ev);
		df4.x -= 1.0 / 2 * third_directional_deriv(gradientField, x + 2, y, z, ev);
	}
	else if (x > 2048 - 2) {
		df4.x += 3.0 / 2 * third_directional_deriv(gradientField, x, y, z, ev);
		df4.x -= 2 * third_directional_deriv(gradientField, x - 1, y, z, ev);
		df4.x += 1.0 / 2 * third_directional_deriv(gradientField, x - 2, y, z, ev);
	}
	else {
		df4.x += third_directional_deriv(gradientField, x + 1, y, z, ev);
		df4.x -= third_directional_deriv(gradientField, x - 1, y, z, ev);
	}

	df4.y = 0;
	if (y < 1) {
		df4.y -= 3.0 / 2 * third_directional_deriv(gradientField, x, y, z, ev);
		df4.y += 2 * third_directional_deriv(gradientField, x, y + 1, z, ev);
		df4.y -= 1.0 / 2 * third_directional_deriv(gradientField, x, y + 2, z, ev);
	}
	else if (y > 2048 - 2) {
		df4.y += 3.0 / 2 * third_directional_deriv(gradientField, x, y, z, ev);
		df4.y -= 2 * third_directional_deriv(gradientField, x, y - 1, z, ev);
		df4.y += 1.0 / 2 * third_directional_deriv(gradientField, x, y - 2, z, ev);
	}
	else {
		df4.y += third_directional_deriv(gradientField, x, y + 1, z, ev);
		df4.y -= third_directional_deriv(gradientField, x, y - 1, z, ev);
	}

	df4.z = 0;
	if (z < 1) {
		df4.z -= 3.0 / 2 * third_directional_deriv(gradientField, x, y, z, ev);
		df4.z += 2 * third_directional_deriv(gradientField, x, y, z + 1, ev);
		df4.z -= 1.0 / 2 * third_directional_deriv(gradientField, x, y, z + 2, ev);
	}
	else if (z > 201 - 2) {
		df4.z += 3.0 / 2 * third_directional_deriv(gradientField, x, y, z, ev);
		df4.z -= 2 * third_directional_deriv(gradientField, x, y, z - 1, ev);
		df4.z += 1.0 / 2 * third_directional_deriv(gradientField, x, y, z - 2, ev);
	}
	else {
		df4.z += third_directional_deriv(gradientField, x, y, z + 1, ev);
		df4.z -= third_directional_deriv(gradientField, x, y, z - 1, ev);
	}

	return df4.x * ev.x + df4.y * ev.y + df4.z * ev.z;
}


void segment_blob(
	std::vector<std::tuple<int, int, int>>& points, 
	std::vector<std::tuple<int,int,int>>& boundary, 
	std::tuple<int, int, int>& local_max, 
	float3* gradientField) {
	int back = 0; 
	if (!points.size() == 0) {
		std::cout << "starting points size is not 0. Returning out of segment_nucleus" << std::endl;
		return;
	}
	char* visited = new char[2048*2048*201];
	memset(visited, 0, sizeof(char) * 2048 * 2048 * 201);
	int x_orig, y_orig, z_orig;	
	std::tie(x_orig, y_orig, z_orig) = local_max;
	//std::cout << "Computing blob at maximum: (" << x_orig << ", " << y_orig << ", " << z_orig << ")" << std::endl; 
	//float max_laplacian = laplacian[x_orig + 2048 * y_orig + 2048 * 2048 * z_orig];
	//if (max_laplacian > 0) { 
		//std::cout << "MAX LAPLACIAN POSITIVE: " << max_laplacian << std::endl;
	//}
	//else {
		//std::cout << "Max Laplacian negative: " << max_laplacian << std::endl;
	//}
	// breadth first search from the local max
	// radial derivative  d^2 f / dr^2 
	points.push_back(std::make_tuple(x_orig, y_orig, z_orig));

	while(points.size() > back) {
		int x, y, z;

		std::tie(x, y, z) = points.at(back);
		back++; 

		std::array<std::tuple<int, int, int>, 6> new_points = {
		std::make_tuple(x + 1, y, z),
		std::make_tuple(x - 1, y, z),
		std::make_tuple(x, y + 1, z),
		std::make_tuple(x, y - 1, z),
		std::make_tuple(x, y, z + 1),
		std::make_tuple(x, y, z - 1)
		};

		// compute radial second derivative.		
		for (int i = 0; i < new_points.size(); i++) {
			// if in points, continue
			int x2, y2, z2;
			std::tie(x2, y2, z2) = new_points.at(i);

			if (x2 - x_orig > 1023 || y2 - y_orig > 1023 || z2 - z_orig > 100 || x2 - x_orig < -1023 || y2 - y_orig < -1023 || z2 - z_orig < -100) {
				std::cout << "Error in segment_blob: cell too big: " << x2-x_orig << " " << y2-y_orig << " " << z2-z_orig << std::endl;
				delete[] visited;
				return;
			}
			
			if (visited[(1024 + x2 - x_orig) + (1024 + y2 - y_orig)*2048 + (100 + z2 - z_orig)*2048*2048] == 1) {
				continue;
			}
			if (x2 == 0 || y2 == 0 || z2 == 0 || x2 == 2047 || y2 == 2047 || z2 == 200) {
				boundary.push_back(std::make_tuple(x2, y2, z2));
			}
			else {

				// instead of using laplacian, we look at whether the rate of change of the radial second derivative is positive.
				// so compute d^2 f / dr^2 at new point, if greater than zero, end.
				// d^2 f / dr^2 = grad ( grad(f) dot r ) dot r

				float sixth_order_centered[4] = { 0, 3.0 / 4, -3.0 / 20, 1.0 / 60 };
				float sixth_order_forward[7] = { -49.0 / 20, 6.0, -15.0 / 2, 20.0 / 3, -15.0 / 4, 6.0 / 5, -1.0 / 6 };
				float gradx = 0, grady = 0, gradz = 0;

				std::array<std::array<float, 3>, 3> hessian = {};
				hessianAt(gradientField, x2, y2, z2, hessian);

				float3 r;
				r.x = x2 - x_orig;
				r.y = y2 - y_orig;
				r.z = z2 - z_orig;
				float r_mag = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
				r.x = r.x / r_mag;
				r.y = r.y / r_mag;
				r.z = r.z / r_mag;

				float3 theta_hat;
				theta_hat.x = -r.y / sqrt(1 - r.z * r.z);
				theta_hat.y = r.x / sqrt(1 - r.z * r.z);
				theta_hat.z = 0;
				float3 phi_hat;
				phi_hat.x = theta_hat.y * r.z;
				phi_hat.y = -theta_hat.x * r.z;
				phi_hat.z = -sqrt(1 - r.z * r.z);


				float3 next_gradient = gradientField[x2 + 2048 * y2 + 2048 * 2048 * z2];
				float next_grad_mag = sqrt(next_gradient.x * next_gradient.x + next_gradient.y * next_gradient.y + next_gradient.z * next_gradient.z);
				float3 grad_dir;
				grad_dir.x = next_gradient.x / next_grad_mag;
				grad_dir.y = next_gradient.y / next_grad_mag;
				grad_dir.z = next_gradient.z / next_grad_mag;

				// what is the second derivative 
				float d2fdg2 = grad_dir.x * (hessian[0][0] * grad_dir.x + hessian[1][0] * grad_dir.y + hessian[2][0] * grad_dir.z) +
					grad_dir.y * (hessian[1][0] * grad_dir.x + hessian[1][1] * grad_dir.y + hessian[2][1] * grad_dir.z) +
					grad_dir.z * (hessian[2][0] * grad_dir.x + hessian[2][1] * grad_dir.y + hessian[2][2] * grad_dir.z);

				// 
				float d2fdtheta2 = theta_hat.x * (hessian[0][0] * theta_hat.x + hessian[1][0] * theta_hat.y + hessian[2][0] * theta_hat.z) +
					theta_hat.y * (hessian[1][0] * theta_hat.x + hessian[1][1] * theta_hat.y + hessian[2][1] * theta_hat.z) +
					theta_hat.z * (hessian[2][0] * theta_hat.x + hessian[2][1] * theta_hat.y + hessian[2][2] * theta_hat.z);

				float d2fdphi2 = phi_hat.x * (hessian[0][0] * phi_hat.x + hessian[1][0] * phi_hat.y + hessian[2][0] * phi_hat.z) +
					phi_hat.y * (hessian[1][0] * phi_hat.x + hessian[1][1] * phi_hat.y + hessian[2][1] * phi_hat.z) +
					phi_hat.z * (hessian[2][0] * phi_hat.x + hessian[2][1] * phi_hat.y + hessian[2][2] * phi_hat.z);

				float d2fdr2 = r.x * (hessian[0][0] * r.x + hessian[1][0] * r.y + hessian[2][0] * r.z) +
					r.y * (hessian[1][0] * r.x + hessian[1][1] * r.y + hessian[2][1] * r.z) +
					r.z * (hessian[2][0] * r.x + hessian[2][1] * r.y + hessian[2][2] * r.z);

				float hessianphitheta[3][2];
				hessianphitheta[0][0] = hessian[0][0] * theta_hat.x + hessian[1][0] * theta_hat.y + hessian[2][0] * theta_hat.z;
				hessianphitheta[1][0] = hessian[1][0] * theta_hat.x + hessian[1][1] * theta_hat.y + hessian[2][1] * theta_hat.z;
				hessianphitheta[2][0] = hessian[2][0] * theta_hat.x + hessian[2][1] * theta_hat.y + hessian[2][2] * theta_hat.z;

				hessianphitheta[0][1] = hessian[0][0] * phi_hat.x + hessian[1][0] * phi_hat.y + hessian[2][0] * phi_hat.z;
				hessianphitheta[1][1] = hessian[1][0] * phi_hat.x + hessian[1][1] * phi_hat.y + hessian[2][1] * phi_hat.z;
				hessianphitheta[2][1] = hessian[2][0] * phi_hat.x + hessian[2][1] * phi_hat.y + hessian[2][2] * phi_hat.z;

				float hessianphitheta2[2][2];
				hessianphitheta2[0][0] = hessianphitheta[0][0] * theta_hat.x + hessianphitheta[1][0] * theta_hat.y + hessianphitheta[2][0] * theta_hat.z;
				hessianphitheta2[0][1] = hessianphitheta[0][1] * theta_hat.x + hessianphitheta[1][1] * theta_hat.y + hessianphitheta[2][1] * theta_hat.z;
				hessianphitheta2[1][0] = hessianphitheta[0][0] * phi_hat.x + hessianphitheta[1][0] * phi_hat.y + hessianphitheta[2][0] * phi_hat.z;
				hessianphitheta2[1][1] = hessianphitheta[0][1] * phi_hat.x + hessianphitheta[1][1] * phi_hat.y + hessianphitheta[2][1] * phi_hat.z;

				float tr = hessianphitheta2[0][0] + hessianphitheta2[1][1];
				float d = hessianphitheta2[0][0] * hessianphitheta2[1][1] - hessianphitheta2[0][1] * hessianphitheta2[1][0];
				float lambda1 = (tr + sqrt(tr * tr - 4 * d)) / 2;
				float lambda2 = (tr - sqrt(tr * tr - 4 * d)) / 2;

				float eig1, eig2, eig3;
				eigenvalues(hessian, eig1, eig2, eig3);
				/*
				// placeholder values if eigs are not positive
				float d4fdeig1_4 = -1, d4fdeig2_4 = -1, d4fdeig3_4 = -1;

				// this gives us the priciple curvature
				float3 ev1, ev2, ev3;
				eigenvector(hessian, eig1, ev1);
				eigenvector(hessian, eig2, ev2);
				eigenvector(hessian, eig3, ev3);

				
				if (eig1 >= 0) { // if eig1 >=0 check if the fourth derivative is >=0 
					d4fdeig1_4 = fourth_directional_deriv(gradientField, x2, y2, z2, ev1);
				}
				if (eig2 >= 0) {
					d4fdeig2_4 = fourth_directional_deriv(gradientField, x2, y2, z2, ev2);
				}
				if (eig3 >= 0) {
					d4fdeig3_4 = fourth_directional_deriv(gradientField, x2, y2, z2, ev3);
				}*/

				if ( //(eig1 >= 0 && d4fdeig1_4 >=0) || (eig2 >= 0 && d4fdeig2_4 >=0) || (eig3 >= 0 && d4fdeig3_4 >= 0)
					(eig1 >= 0) || (eig2 >= 0) || (eig3 >= 0 )
					) {
					//|| d2fdr2 >= 0 || lambda1 >= 0 || lambda2 >= 0
					//d2fdg2 >=0
					//((eig1 >=0) &&  (eig2>= 0)) || ((eig1>=0) && (eig3>=0)) || ((eig2>=0)&&(eig3>=0)) 
					//eig1+eig2+eig3 >= 0
					//|| d2fdr2 >= 0 
					//gaussian[x2 + 2048 * y2 + 2048 * 2048 * z2] < 1)
					//|| next_laplacian - d2fdr2 >= 0 // || next_laplacian < this_laplacian //-1/(sigma*sigma)*0.37
					//|| next_gradient.x * r.x + next_gradient.y * r.y + next_gradient.z * r.z >= 0
					//|| next_gaussian > this_gaussian // solves kissing problem?

					// if laplacian is greater than or equal 0, then we're in the boundary.
					// also if on boundary of image
					boundary.push_back(std::make_tuple(x2, y2, z2));
				}
				else {
					// if laplacian is less than 0, then we're in points
					points.push_back(std::make_tuple(x2, y2, z2));
				}
			}
			visited[(1024 + x2 - x_orig) + (1024 + y2 - y_orig) * 2048 + (100 + z2 - z_orig) * 2048 * 2048] = 1;
		}
	}
	delete[] visited;
}


void medianFilter3x3(uint16_t* voxels, uint16_t* filtered) {
	concurrency::parallel_for(0, 201, [&voxels, &filtered](int z) {
		for (int y = 0; y < 2048; y++) {
			for (int x = 1; x < 2048; x++) {
				if (z == 0 || z == 200 || y == 0 || y == 2047 || x == 0 || x == 2047) {
					filtered[x + 2048 * y + 2048 * 2048 * z] = 0;
					continue;
				}
				uint16_t A[27], B[27];
				uint16_t* current;
				uint16_t* next;
				current = A;
				next = B;
				for (int i = -1; i < 2; i++) {
					for (int j = -1; j < 2; j++) {
						for (int k = -1; k < 2; k++) {
							current[(i + 1) + 3 * (j + 1) + 3 * 3 * (k + 1)] = voxels[(x + k) + 2048 * (y + j) + 2048 * 2048 * (z + i)];
						}
					}
				}
				// find nth member of the list
				int start_index = 0;
				int len = 27;
				int n = 14;
				while (len > 0) {
					int num_lesser = 0, num_greater = 0;
					uint16_t pivot = current[start_index];
					for (int i = 1; i < len; i++) {
						if (current[start_index + i] > pivot) {
							next[start_index + len - 1 - num_greater] = current[start_index + i];
							num_greater++;
						}
						else {
							next[start_index + num_lesser] = current[start_index + i];
							num_lesser++;
						}
					}
					if (num_lesser == n - 1) { // we've found the nth member
						filtered[x + 2048 * y + 2048 * 2048 * z] = pivot;
						len = 0;
					}
					else if (num_lesser > n - 1) { // n-1 is in the lesser list
						next[start_index + len - 1 - num_greater] = pivot;
						num_greater++;
						len = len - num_greater;
					}
					else { // n-1 is in the greater list
						next[start_index + num_lesser] = pivot;
						num_lesser++;
						start_index += num_lesser;
						len = len - num_lesser;
						n = n - num_lesser;
					}
					uint16_t* temp = current;
					current = next;
					next = temp;
				}
			}
		}
		});
	
}

 // medium fast gaussian filter
 // In this function, for optimization purposes I try to use some pointer
 // arithmetic so that I'm always reading from continuous memory.
void gaussian_filter3D_parallel(uint16_t* input, int sigmaxy, int sigmaz, float* result) {
	float float_sigmaxy = (float)sigmaxy;
	float float_sigmaz = (float)sigmaz;
	// prerun expensive exp operation in 1d array
	float* kernelxy = new float[2 * GAUSSIAN_CUTOFF * sigmaxy];
	float* kernelz = new float[2 * GAUSSIAN_CUTOFF * sigmaz];
	float normxy = 1.0 / sqrt(2 * M_PI * float_sigmaxy * float_sigmaxy);
	float normz = 1.0 / sqrt(2 * M_PI * float_sigmaz * float_sigmaz);

	for (int r = 0; r <= 2*GAUSSIAN_CUTOFF*sigmaxy; r++) {
		kernelxy[r] = exp(-r * r / (2 * float_sigmaxy * float_sigmaxy)) * normxy;
	}
	for (int r = 0; r <= 2*GAUSSIAN_CUTOFF*sigmaz; r++) {
		kernelz[r] = exp(-r * r / (2 * float_sigmaz * float_sigmaz)) * normz;
	}

	float* temp = new float[2048 * 2048 * 201];

	concurrency::static_partitioner partitioner;
	//concurrency::critical_section cs;

	concurrency::parallel_for(0, 2048, [&input, &result, &sigmaxy, &kernelxy](int x) {
	//for (int x = 0; x < 2048; x++) {
		for (int z = 0; z < 201; z++) {
			// copy to contiguous memory for better speed.
			uint16_t yaxis[2048];
			for (int y = 0; y < 2048; y++) {
				yaxis[y] = input[x + 2048 * 2048 * z + y * 2048];
			}

			for (int y = 0; y < 2048; y++) {
				float sum = 0;
				int j_min, j_max;
				if (y - GAUSSIAN_CUTOFF * sigmaxy > 0) {
					j_min = y - GAUSSIAN_CUTOFF * sigmaxy;
				}
				else {
					j_min = 0;
				}
				if (y + GAUSSIAN_CUTOFF * sigmaxy < 2048) {
					j_max = y + GAUSSIAN_CUTOFF * sigmaxy;
				}
				else {
					j_max = 2048;
				}

				for (int j = j_min; j < y; j++) {
					sum += yaxis[j] * kernelxy[y - j];
				}
				for (int j = y; j < j_max; j++) {
					sum += yaxis[j] * kernelxy[j - y];
				}
				result[z * 2048 + 2048 * 201 * x + y] = sum;
			}
		}
		}, partitioner);
	//}

	concurrency::parallel_for(0, 2048, [&result, &temp, &sigmaz, &kernelz](int y) {
			//for (int y = 0; y < 2048; y++) {
				for (int x = 0; x < 2048; x++) {
					// copy to contiguous memory
					float zaxis[201];
					for (int z = 0; z < 201;z++) {
						zaxis[z] = result[y + 2048 * 201 * x + z * 2048];
					}
					for (int z = 0; z < 201; z++) {
						int j_min, j_max;
						float sum = 0;
						if (z - GAUSSIAN_CUTOFF * sigmaz > 0) {
							j_min = z - GAUSSIAN_CUTOFF * sigmaz;
						}
						else {
							j_min = 0;
						}
						if (z + GAUSSIAN_CUTOFF * sigmaz < 201) {
							j_max = z + GAUSSIAN_CUTOFF * sigmaz;
						}
						else {
							j_max = 201;
						}

						for (int j = j_min; j < z; j++) {
							sum += zaxis[j] * kernelz[z - j];
						}
						for (int j = z; j < j_max; j++) {
							sum += zaxis[j] * kernelz[j - z];
						}
						temp[201 * x + 2048 * 201 * y + z] = sum;
					}
				}
				}, partitioner);
			//}

			concurrency::parallel_for(0, 201, [&temp, &result, &sigmaxy, &kernelxy](int z) {
			//for (int z = 0; z < 201; z++) {
				for (int y = 0; y < 2048; y++) {
					int j_min, j_max;
					float xaxis[2048];
					for (int x = 0; x < 2048; x++) {
						xaxis[x] = temp[z + 2048 * 201 * y + x * 201];
					}
						for (int x = 0; x < 2048; x++) {
							float sum = 0;
							if (x - GAUSSIAN_CUTOFF * sigmaxy > 0) {
								j_min = x - GAUSSIAN_CUTOFF * sigmaxy;
							}
							else {
								j_min = 0;
							}
							if (x + GAUSSIAN_CUTOFF * sigmaxy < 2048) {
								j_max = x + GAUSSIAN_CUTOFF * sigmaxy;
							}
							else {
								j_max = 2048;
							}
							for (int j = j_min; j < x; j++) {
								sum += xaxis[j] * kernelxy[x - j];
							}
							for (int j = x; j < j_max; j++) {
								sum += xaxis[j] * kernelxy[j - x];
							}
							//cs.lock();
							result[2048 * y + 2048 * 2048 * z + x] = sum;
							//cs.unlock();
						}
					}
						}, partitioner);

	// clean up
	// delete 
	delete[] temp;
	delete[] kernelxy;
	delete[] kernelz;
}


// Loading tiff data into a packed array form
// at this point not flexible to dimensions of array.
void load_tiff2(const char* filename, uint16_t* output) {
	TIFF* tiff = TIFFOpen(filename, "r");

	// turn off warnings
	TIFFSetWarningHandler(0);

	int zed = 0;
	do {
		int ied = 0;
		for (int i = 0; i < 2048; i++) {
			TIFFReadScanline(tiff, (output + ied + zed), i);
			ied += 2048;
		}
		zed += 2048 * 2048;
	} while (TIFFReadDirectory(tiff));
	TIFFClose(tiff);
}

int main()
{
	// files
	std::array<char*, 4> construct_dirs = {
		"MouseBrain2021-10-21\\0x\\0x_slice1_",
		"MouseBrain2021-10-21\\1x\\1x_",
		"MouseBrain2021-10-21\\10x\\10x_",
		"MouseBrain2021-10-21\\endogenous_0x\\endogenous_0x_"
	};
	//FILE* dots_csv;
	//FILE* nucleii_csv;
	//
	FILE* dots_csv = fopen("dots.csv", "w");
	fprintf(dots_csv, "id,x,y,z,cn_id,cn_x,cn_y,cn_z,channel,dir,pos\n");
	fclose(dots_csv);

	FILE* nucleii_csv = fopen("nucleii.csv", "w");
	fprintf(nucleii_csv, "id,x,y,z,size,validSize,sum488,dir,pos\n");
	fclose(nucleii_csv);

	for (char* dir : construct_dirs) {
		int pos = 1;
		
		char dapi_filename[500];
		sprintf(dapi_filename, "%s#%d_New.tif", dir, pos);

		FILE* dapi_file;
		while ((dapi_file = fopen(dapi_filename, "r")) && pos < 10) {
			fclose(dapi_file);
			std::cout << "Working on " << dapi_filename << std::endl;

			// fill out filenames
			char file_488[500];
			char file_594[500];
			char file_647[500];
			sprintf(file_488, "%s#%d_New_1.tif", dir, pos);
			sprintf(file_594, "%s#%d_New_2.tif", dir, pos);
			sprintf(file_647, "%s#%d_New_4.tif", dir, pos);


			int threshold_405_lower = 120;
			int threshold_405_higher = 171;
			//std::cout << "please observe this dog..." << std::endl;
			std::cout << "Loading DAPI...";
			//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
			uint16_t* stack = new uint16_t[2048 * 2048 * 201];
			load_tiff2(dapi_filename, stack);
			std::cout << " done." << std::endl;
			time_t start, end;

			uint16_t* filtered = new uint16_t[2048 * 2048 * 201];
			std::cout << "Computing median filter...";
			time(&start);
			medianFilter3x3(stack, filtered);
			time(&end);
			std::cout << "done." << std::endl;
			std::cout << "median filter took " << end - start << "seconds" << std::endl;

			uint16_t* pointer = filtered;
			for (int z = 0; z < 201; z++) {
				for (int y = 0; y < 2048; y++) {
					for (int x = 0; x < 2048; x++) {
						if ((*pointer) <= threshold_405_lower) {
							*pointer = 0;
						}
						else if ((*pointer) >= threshold_405_higher) {
							*pointer = threshold_405_higher - threshold_405_lower;
						}
						else {
							*pointer = (*pointer) - threshold_405_lower;
						}
						pointer++;
					}
				}
			}

			//cv::waitKey(0);//
			// Gaussian filter??
			//cv::Mat* new_mat = new cv::Mat(mats->at(0)->size().height, mats->at(0)->size().width, CV_32F);
			std::cout << "Computing gaussian filter...";
			time(&start);
			int sigma = 20;
			float* gaussian_result = new float[2048 * 2048 * 201];
			gaussian_filter3D_parallel(filtered, sigma, 10, gaussian_result);
			std::cout << "done." << std::endl;

			time(&end);
			std::cout << "gaussian filter took " << difftime(end, start) << " seconds." << std::endl;

			std::vector<std::tuple<int, int, int>> maxima;
			time(&start);
			findMaxima(gaussian_result, maxima);
			time(&end);
			std::cout << "Found " << maxima.size() << " nuclei." << std::endl;
			std::cout << "findMaxima took " << difftime(end, start) << " seconds." << std::endl;

			//std::cout << "Redoing gaussian filter with larger sigma to eliminate heterochromatic puncta";
			//gaussian_result = gaussian_filter3D_parallel(filtered, sigma*2, sigma);
			//std::cout << "done." << std::endl;

			//float* laplacian = new float[2048 * 2048 * 201];
			//std::cout << "computing laplacian...";             
			//time(&start);
			//laplacianFilter3D(gaussian_result, laplacian);
			//std::cout << "done." << std::endl; 
			//time(&end);
			//std::cout << "Laplacian filter took: " << difftime(end, start) << " seconds." << std::endl;

			int gradientfieldsize = 2048 * 2048 * 201;
			float3* gradientField = new float3[gradientfieldsize];
			std::cout << "computing gradient field...";
			time(&start);
			gradientField3d(gaussian_result, gradientField);
			std::cout << "done." << std::endl;
			time(&end);
			std::cout << "Gradient field took: " << difftime(end, start) << " seconds." << std::endl;
			delete[] gaussian_result;

			//std::cout << "Testing gradient field" << std::endl;
			//std::cout << "test test: " << gradientField[2047 + 2048 * 2047 + 2048 * 2048 * 200].y << std::endl;
			//return 0;

			std::cout << "segmenting blobs...";
			std::vector<Nucleus*> nucleii;
			nucleii.resize(maxima.size());
			time(&start);
			concurrency::static_partitioner partitioner;

			concurrency::parallel_for(0, (int)maxima.size(), [&maxima, &nucleii, gradientField](int j) {
				Nucleus* blob = new Nucleus;
				blob->id = j;
				blob->points.reserve(200000);
				blob->boundary.reserve(20000);
				blob->local_max = maxima.at(j);
				segment_blob(blob->points, blob->boundary, maxima.at(j), gradientField);
				nucleii.at(j) = blob;
				//std::cout << "blob has " << blob->points.size() << " voxels" << std::endl;
				//std::cout << "boundary has " << blob->boundary.size() << " voxels" << std::endl;
				}, partitioner);
			time(&end);

			std::cout << "done." << std::endl;
			std::cout << "blob segmentation took: " << difftime(end, start) << " seconds." << std::endl;

			int num_good = 0;
			for (int i = 0; i < nucleii.size(); i++) {
				if (nucleii.at(i)->validSize()) {
					num_good++;
				}
			}
			std::cout << "Found " << num_good << " out of " << nucleii.size() << " maxima to be good..." << std::endl;

			std::cout << "Loading 488...";
			uint16_t* stack488 = new uint16_t[2048 * 2048 * 201];
			load_tiff2(file_488, stack488);
			std::cout << "done." << std::endl;
			/*
			pointer = stack488;
			for (int z = 0; z < 201; z++) {
				for (int y = 0; y < 2048; y++) {
					for (int x = 0; x < 2048; x++) {
						if ((*pointer) <= 160) {
							*pointer = 0;
						}
						else {
							*pointer = (*pointer) - 160;
						}
						pointer++;
					}
				}
			}*/

			//std::cout << "median filter 488...";
			//uint16_t* median488 = new uint16_t[2048 * 2048 * 201];
			//medianFilter3x3(stack488, median488);
			//std::cout << "done." << std::endl;

			//std::cout << "please observe this dog..." << std::endl;
			std::cout << "Loading 594...";
			//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
			uint16_t* stack594 = new uint16_t[2048 * 2048 * 201];
			load_tiff2(file_594, stack594);
			std::cout << " done." << std::endl;

			std::cout << "median 594...";
			uint16_t* median594 = new uint16_t[2048 * 2048 * 201];
			medianFilter3x3(stack594, median594);
			std::cout << "done." << std::endl;
			delete[] stack594;

			pointer = median594;
			for (int z = 0; z < 201; z++) {
				for (int y = 0; y < 2048; y++) {
					for (int x = 0; x < 2048; x++) {
						if ((*pointer) <= 135) {
							*pointer = 0;
						}
						else {
							*pointer = (*pointer) - 135;
						}
						pointer++;
					}
				}
			}
			std::cout << "computing 594 gaussian...";
			float* gaussian_594 = new float[2048 * 2048 * 201];
			gaussian_filter3D_parallel(median594, 2, 1, gaussian_594);
			std::cout << "done." << std::endl;

			std::vector<std::tuple<int, int, int>> maxima594;
			findMaxima(gaussian_594, maxima594);
			std::cout << "Found " << maxima594.size() << " dots." << std::endl;

			std::cout << "computing 594 gradient field...";
			time(&start);
			gradientField3d(gaussian_594, gradientField);
			std::cout << "done." << std::endl;
			time(&end);
			std::cout << "Gradient field took: " << difftime(end, start) << " seconds." << std::endl;
			delete[] gaussian_594;

			std::cout << "segmenting 594 dots...";
			std::vector<Dot*> dots594;
			dots594.resize(maxima594.size());
			concurrency::parallel_for(0, (int)maxima594.size(), [&maxima594, &gradientField, &dots594](int j) {
				Dot* blob = new Dot;
				blob->id = j;
				blob->points.reserve(2000);
				blob->boundary.reserve(200);
				blob->local_max = maxima594.at(j);
				segment_blob(blob->points, blob->boundary, maxima594.at(j), gradientField);
				dots594.at(j) = blob;
				}, partitioner);
			std::cout << "done." << std::endl;
			std::cout << "Segmented " << dots594.size() << " dots." << std::endl;

			//std::cout << "please observe this dog..." << std::endl;
			std::cout << "Loading 640...";
			//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
			uint16_t* stack647 = new uint16_t[2048 * 2048 * 201];
			load_tiff2(file_647, stack647);
			std::cout << " done." << std::endl;

			std::cout << "median 640...";
			uint16_t* median640 = new uint16_t[2048 * 2048 * 201];
			medianFilter3x3(stack647, median640);
			std::cout << "done." << std::endl;
			delete[] stack647;

			// analyzing pos2 on 0x seems to indicate that 220 is the correct number here.
			pointer = median640;
			for (int z = 0; z < 201; z++) {
				for (int y = 0; y < 2048; y++) {
					for (int x = 0; x < 2048; x++) {
						if ((*pointer) <= 500) {
							*pointer = 0;
						}
						else {
							*pointer = (*pointer) - 500;
						}
						pointer++;
					}
				}
			}
			std::cout << "computing 647 gaussian...";
			float* gaussian_647 = new float[2048 * 2048 * 201];
			gaussian_filter3D_parallel(median640, 2, 1, gaussian_647);
			std::cout << "done." << std::endl;

			std::vector<std::tuple<int, int, int>> maxima647;
			findMaxima(gaussian_647, maxima647);
			std::cout << "Found " << maxima647.size() << " dots." << std::endl;

			std::cout << "computing 647 gradient field...";
			time(&start);
			gradientField3d(gaussian_647, gradientField);
			std::cout << "done." << std::endl;
			time(&end);
			std::cout << "Gradient field took: " << difftime(end, start) << " seconds." << std::endl;
			delete[] gaussian_647;

			std::cout << "segmenting 640 dots...";
			std::vector<Dot*> dots647;
			dots647.resize(maxima647.size());
			concurrency::parallel_for(0, (int)maxima647.size(), [&maxima647, &gradientField, &dots647](int j) {
				Dot* blob = new Dot;
				blob->id = j;
				blob->points.reserve(2000);
				blob->boundary.reserve(200);
				blob->local_max = maxima647.at(j);
				segment_blob(blob->points, blob->boundary, maxima647.at(j), gradientField);
				dots647.at(j) = blob;
				}, partitioner);
			std::cout << "done." << std::endl;
			
			
			nucleii_csv = fopen("nucleii.csv", "a");
			for (int i = 0; i < maxima.size(); i++) {
				Nucleus* nuc = nucleii.at(i);
				int max_x, max_y, max_z;
				std::tie(max_x, max_y, max_z) = nuc->local_max;
				int sumgreen = 0;
				for (std::tuple<int, int, int> pos : nuc->points) {
					int xg, yg, zg;
					std::tie(xg, yg, zg) = pos; 
					sumgreen += stack488[xg + 2048 * yg + 2048 * 2048 * zg];
				}
				fprintf(nucleii_csv, "%d,%d,%d,%d,%d,%s,%d,%s,%d\n", nuc->id, max_x, max_y, max_z, nuc->points.size(), nuc->validSize() ? "TRUE" : "FALSE", sumgreen, dir, pos);
			}
			fclose(nucleii_csv);
			
			// make the dots csv
			// obviously we need x,y,z,channel, 
			// also closest nuclei id, and x,y,z
			dots_csv = fopen("dots.csv", "a");
			for (int i = 0; i < dots594.size(); i++) {
				Dot* dot = dots594.at(i);
				if (!dot->validSize()) continue;
				int x, y, z;
				std::tie(x, y, z) = dot->local_max;

				int cn_x, cn_y, cn_z, cn_id = 0;
				Nucleus* nuc = nucleii.at(0);
				std::tie(cn_x, cn_y, cn_z) = nuc->local_max;

				float closest = (x - cn_x) * (x - cn_x) + (y - cn_y) * (y - cn_y) + (z - cn_z) * (z - cn_z);
				for (int j = 1; j < nucleii.size(); j++) {
					if (!nucleii.at(j)->validSize()) continue;

					int n_x, n_y, n_z;
					std::tie(n_x, n_y, n_z) = nucleii.at(j)->local_max;
					float this_dist = (x - n_x) * (x - n_x) + (y - n_y) * (y - n_y) + (z - n_z) * (z - n_z);
					if (this_dist < closest) {
						nuc = nucleii.at(j);
						closest = this_dist;
						cn_x = n_x;
						cn_y = n_y;
						cn_z = n_z;
						cn_id = j;
					}
				}
				if (closest < 125 * 125) {
					nuc->close_dots594.push_back(dot);
					fprintf(dots_csv, "%d,%d,%d,%d,%d,%d,%d,%d,594,%s,%d\n", i, x, y, z, cn_id, cn_x, cn_y, cn_z, dir, pos);
				}
			}

			for (int i = 0; i < dots647.size(); i++) {
				Dot* dot = dots647.at(i);
				int x, y, z;
				std::tie(x, y, z) = dot->local_max;

				int cn_x, cn_y, cn_z, cn_id = dots594.size();
				Nucleus* nuc = nucleii.at(0);
				std::tie(cn_x, cn_y, cn_z) = nuc->local_max;

				float closest = (x - cn_x) * (x - cn_x) + (y - cn_y) * (y - cn_y) + (z - cn_z) * (z - cn_z);
				for (int j = 1; j < nucleii.size(); j++) {
					if (!nucleii.at(j)->validSize()) continue;
					int n_x, n_y, n_z;
					std::tie(n_x, n_y, n_z) = nucleii.at(j)->local_max;
					float this_dist = (x - n_x) * (x - n_x) + (y - n_y) * (y - n_y) + (z - n_z) * (z - n_z);
					if (this_dist < closest) {
						nuc = nucleii.at(j);
						closest = this_dist;
						cn_x = n_x;
						cn_y = n_y;
						cn_z = n_z;
						cn_id = j;
					}
				}
				if (closest < 125 * 125) {
					nuc->close_dots640.push_back(dot);
					fprintf(dots_csv, "%d,%d,%d,%d,%d,%d,%d,%d,640,%s,%d\n", i + maxima594.size(), x, y, z, cn_id, cn_x, cn_y, cn_z, dir, pos);
				}
			}
			fclose(dots_csv);
			std::cout << "Wrote dots to dots.csv" << std::endl;

			/*
			int i = 0;
			int k = 0;
			bool blink = true;
			int key = 0;
			do {
				//cv::Mat img(2048, 2048, CV_32F, gaussian_result + 2048*2048*i);
				cv::Mat img(2048, 2048, CV_16U, stack + 2048 * 2048 * i);
				cv::Mat dst(2048, 2048, CV_16UC3);
				cv::cvtColor(img, dst, cv::COLOR_GRAY2BGR);

				for (int l = 0; l < nucleii.size(); l++) {
					Nucleus* nuc = nucleii.at(l);
					std::tuple<int, int, int> max = nuc->local_max;

					if (nuc->validSize()) {
						int z = std::get<2>(max);
						if (abs(i - z) < 10) {
							int cx = std::get<0>(max);
							int cy = std::get<1>(max);
							cv::rectangle(dst, cv::Point(cx - 5, cy - 5), cv::Point(cx + 5, cy + 5), cv::Scalar(0, 0, 65535), cv::FILLED, cv::LINE_8);
						}

						for (int j = 0; j < nuc->boundary.size(); j++) {
							int x, y, z;

							std::tie(x, y, z) = nuc->boundary.at(j);

							if (z == i) {
								if (l == k) {
									BGR& bgr = dst.ptr<BGR>(y)[x];
									bgr.red = 255;
									bgr.green = 0;
									bgr.blue = 0;
								}
								else {
									BGR& bgr = dst.ptr<BGR>(y)[x];
									bgr.red = 0;
									bgr.green = 255;
									bgr.blue = 0;
								}
							}
						}
					}
				}
				//cv::Mat resized(512, 512, CV_32F);
				cv::Mat resized(512, 512, CV_16UC3);
				cv::resize(dst, resized, cv::Size(512, 512));
				//cv::imshow("Display window", resized);
				cv::imshow("Display window", resized * 255);
				key = cv::waitKey(0);
				if (key == 's') {
					i = (i + 1) % 201;
				}
				else if (key == 'w') {
					i = (201 + i - 1) % 201;
				}
				else if (key == 'a') {
					k = (k + 1) % nucleii.size();
					i = std::get<2>(nucleii.at(k)->local_max);
				}
				else if (key == 'd') {
					k = (nucleii.size() + k - 1) % nucleii.size();
					i = std::get<2>(nucleii.at(k)->local_max);
				}
				else {
					//blink = !blink;
				}

			} while (key != 27);*/

			/*
			i = 0;
			k = 0;
			do {
				cv::Mat img(2048, 2048, CV_32F, gaussian_result + 2048 * 2048 * i);
				double min, max;
				cv::minMaxLoc(img, &min, &max);
				img = (img - min) / (max - min);

				cv::Mat dst(2048, 2048, CV_32FC3);
				cv::cvtColor(img, dst, cv::COLOR_GRAY2BGR);

				//cv::Mat img(2048, 2048, CV_16U, stack + 2048 * 2048 * i);
				for (int l = 0; l < nucleii.size(); l++) {
					Nucleus* nuc = nucleii.at(l);
					std::tuple<int, int, int> max = nuc->local_max;

					if (nuc->validSize()) {
						int z = std::get<2>(max);
						if (abs(i - z) < 10) {
							int cx = std::get<0>(max);
							int cy = std::get<1>(max);
							cv::rectangle(dst, cv::Point(cx - 5, cy - 5), cv::Point(cx + 5, cy + 5), cv::Scalar(0, 0, 255), cv::FILLED, cv::LINE_8);
						}

						for (int j = 0; j < nuc->boundary.size(); j++) {
							int x, y, z;

							std::tie(x, y, z) = nuc->boundary.at(j);

							if (z == i) {
								if (l == k) {
									BGR_float& bgr = dst.ptr<BGR_float>(y)[x];
									bgr.red = 1;
									bgr.green = 0;
									bgr.blue = 0;
								}
								else {
									BGR_float& bgr = dst.ptr<BGR_float>(y)[x];
									bgr.red = 0;
									bgr.green = 1;
									bgr.blue = 0;
								}
							}
						}
					}
				}

				cv::Mat resized(512, 512, CV_32F);
				//cv::Mat resized(512, 512, CV_16U);
				cv::resize(dst, resized, cv::Size(512, 512));

				//cv::imshow("Display window", resized / 20);
				cv::imshow("Display window", resized);
				key = cv::waitKey(0);
				if (key == 's') {
					i = (i + 1) % 201;
				}
				else if (key == 'w') {
					i = (201 + i - 1) % 201;
				}
				else if (key == 'a') {
					k = (k + 1) % nucleii.size();
					i = std::get<2>(nucleii.at(k)->local_max);
				}
				else if (key == 'd') {
					k = (nucleii.size() + k - 1) % nucleii.size();
					i = std::get<2>(nucleii.at(k)->local_max);
				}
				else {
					//blink = !blink;
				}

			} while (key != 27);*/
			

			// on this screen we want to display a single cell with its dots
			
			const int blue_pixel_lower = 120;
			const int blue_pixel_upper = 200;

			const int green_pixel_lower = 0;
			const int green_pixel_upper = 1;

			const int red_pixel_lower = 0;
			const int red_pixel_upper = 1;

			const int white_pixel_lower = 0;
			const int white_pixel_upper = 1;
			/*
			k = 0;
			i = 0;
			do {
				int nx, ny, nz;
				Nucleus* nuc = nucleii.at(k);
				std::tie(nx, ny, nz) = nuc->local_max;
				cv::Mat img(401, 401, CV_16UC3, cv::Scalar(0, 0, 0));

				// we are going to look at range from maxima-150 to maxima+150
				for (int rx = -200; rx <= 200; rx++) {
					if (nx + rx < 0 || nx + rx > 2047) {
						continue;
					}
					for (int ry = -200; ry <= 200; ry++) {
						if (ny + ry < 0 || ny + ry > 2047) {
							continue;
						}
						BGR& bgr = img.ptr<BGR>(ry + 200)[rx + 200];

						uint16_t green = median594[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * i];
						if (green < green_pixel_lower) bgr.green = 0;
						else if (green > green_pixel_upper) bgr.green = 255 * 255;
						else bgr.green = (green - green_pixel_lower) * (255 * 255 / green_pixel_upper);

						uint16_t blue = stack[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * i];
						if (blue < blue_pixel_lower) bgr.blue = 0;
						else if (blue > blue_pixel_upper) bgr.blue = 255 * 255;
						else bgr.blue = (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper);

						uint16_t red = median640[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * i];
						if (red < red_pixel_lower) bgr.red = 0;
						else if (red > red_pixel_upper) bgr.red = 255 * 255;
						else bgr.red = (red - red_pixel_lower) * (255 * 255 / red_pixel_upper);


						/*uint16_t white = stack488[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * i];
						if (white <= white_pixel_lower) {
							//bgr.red = 
						}
						else if (white > white_pixel_upper) {
							bgr.red = 255 * 255;
							bgr.green = 255 * 255;
							bgr.blue = 255 * 255;
						}
						else {
							bgr.red = (white - white_pixel_lower) * (255 * 255 / white_pixel_upper);
							bgr.green = bgr.red;
							bgr.blue = bgr.red;
						} */
					//}
				//}
				//cv::Mat channels[3];
				//cv::split(img, channels);

				//cv::Mat greenMat, redMat;
				//channels[1].convertTo(greenMat, CV_32F);
				//channels[2].convertTo(redMat, CV_32F);

				//cv::Point2d point = cv::phaseCorrelate(greenMat, redMat);
				//std::cout << "phase correlate: (" << point.x << ", " << point.y << ")" << std::endl;
				/*
				for (int j = 0; j < nuc->boundary.size(); j++) {
					int x, y, z;

					std::tie(x, y, z) = nuc->boundary.at(j);

					if (x - nx > 200 || x - nx < -200 || y - ny > 200 || y - ny < -200) continue;
					if (z == i) {
						BGR& bgr = img.ptr<BGR>(y - ny + 200)[x - nx + 200];
						bgr.red = 255 * 255;
						bgr.green = 255 * 255;
						bgr.blue = 255 * 255;

					}
				}

				// for all the dots, show their outline
				for (int j = 0; j < nuc->close_dots594.size(); j++) {
					int dx, dy, dz;
					std::tie(dx, dy, dz) = nuc->close_dots594.at(j)->local_max;

					if (((dz - i) < 5 && i - dz < 5)) {
						// if we are close enough 
						cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(0, 255 * 255, 0), 1, 8, 0);
						cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);

					}

				}

				// for all the dots, show their outline
				for (int j = 0; j < nuc->close_dots640.size(); j++) {
					int dx, dy, dz;
					std::tie(dx, dy, dz) = nuc->close_dots640.at(j)->local_max;

					if (((dz - i) < 5 && i - dz < 5)) {
						//cv::rectangle(img, cv::Point(dx - 5, dy - 5), cv::Point(dx + 5, dy + 5), cv::Scalar(0, 0, 255), cv::FILLED, cv::LINE_8);
						cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(0, 0, 255 * 255), 1, 8, 0);
						cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);
					}

				}
				cv::imshow("Display window", img);
				key = cv::waitKey(0);
				if (key == 's') {
					i = (i + 1) % 201;
				}
				else if (key == 'w') {
					i = (201 + i - 1) % 201;
				}
				else if (key == 'a') {
					do {
						k = (k + 1) % nucleii.size();
					} while (!nucleii.at(k)->validSize());
					i = std::get<2>(nucleii.at(k)->local_max);
				}
				else if (key == 'd') {
					do {
						k = (nucleii.size() + k - 1) % nucleii.size();
					} while (!nucleii.at(k)->validSize());
					i = std::get<2>(nucleii.at(k)->local_max);
				}
				else {
					//blink = !blink;
				}
			} while (key != 27);*/

			// 2d maxproject
			/*k = 0;
			do {
				int nx, ny, nz;
				Nucleus* nuc = nucleii.at(k);
				std::tie(nx, ny, nz) = nuc->local_max;
				i = nz;
				cv::Mat img(401, 401, CV_16UC3, cv::Scalar(0, 0, 0));

				// we are going to look at range from maxima-150 to maxima+150
				for (int rx = -200; rx <= 200; rx++) {
					if (nx + rx < 0 || nx + rx > 2047) {
						continue;
					}
					for (int ry = -200; ry <= 200; ry++) {
						if (ny + ry < 0 || ny + ry > 2047) {
							continue;
						}
						BGR& bgr = img.ptr<BGR>(ry + 200)[rx + 200];

						for (int zh = i - 50; zh < i + 50; zh++) {
							if (zh < 0) continue;
							if (zh > 200) break;
							uint16_t green = median594[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * zh];
							if (bgr.green < 255 * 255 && bgr.green < (green - green_pixel_lower) * (255 * 255 / green_pixel_upper)) {
								if (green > green_pixel_upper) bgr.green = 255 * 255;
								else bgr.green = (green - green_pixel_lower) * (255 * 255 / green_pixel_upper);
							}
							uint16_t blue = stack[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * zh];

							if (bgr.blue < 255 * 255 && bgr.blue < (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper)) {
								if (blue > blue_pixel_upper) bgr.blue = 255 * 255;
								else bgr.blue = (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper);
							}

							uint16_t red = median640[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * zh];
							if (bgr.red < 255 * 255 && bgr.red < (red - red_pixel_lower) * (255 * 255 / red_pixel_upper)) {
								if (red > red_pixel_upper) bgr.red = 255 * 255;
								else bgr.red = (red - red_pixel_lower) * (255 * 255 / red_pixel_upper);
							}
						}
					}
				}

				// for all the dots, show their outline
				for (int j = 0; j < nuc->close_dots594.size(); j++) {
					int dx, dy, dz;
					std::tie(dx, dy, dz) = nuc->close_dots594.at(j)->local_max;
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(0, 255 * 255, 0), 1, 8, 0);
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);

				}

				// for all the dots, show their outline
				for (int j = 0; j < nuc->close_dots640.size(); j++) {
					int dx, dy, dz;
					std::tie(dx, dy, dz) = nuc->close_dots640.at(j)->local_max;
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(0, 0, 255 * 255), 1, 8, 0);
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);
				}

				for (int j = 0; j < nuc->boundary.size(); j++) {
					int x, y, z;

					std::tie(x, y, z) = nuc->boundary.at(j);

					if (x - nx > 200 || x - nx < -200 || y - ny > 200 || y - ny < -200) continue;
					if (z == i) {
						BGR& bgr = img.ptr<BGR>(y - ny + 200)[x - nx + 200];
						bgr.red = 255 * 255;
						bgr.green = 255 * 255;
						bgr.blue = 255 * 255;

					}
				}

				cv::imshow("Display window", img);

				key = cv::waitKey(0);
				if (key == 'a') {
					do {
						k = (k + 1) % nucleii.size();
					} while (!nucleii.at(k)->validSize());
				}
				else if (key == 'd') {
					do {
						k = (nucleii.size() + k - 1) % nucleii.size();
					} while (!nucleii.at(k)->validSize());
				}

			} while (key != 27); */


			// save 2D maxprojects
			// 2d maxproject 

			std::cout << "writing maxprojects...";
			for (int k = 0; k < nucleii.size(); k++) {
				int nx, ny, nz;
				Nucleus* nuc = nucleii.at(k);
				if (!nuc->validSize()) continue;
				std::tie(nx, ny, nz) = nuc->local_max;
				int i = nz;
				cv::Mat img(401, 401, CV_16UC3, cv::Scalar(0, 0, 0));

				// we are going to look at range from maxima-150 to maxima+150
				for (int rx = -200; rx <= 200; rx++) {
					if (nx + rx < 0 || nx + rx > 2047) {
						continue;
					}
					for (int ry = -200; ry <= 200; ry++) {
						if (ny + ry < 0 || ny + ry > 2047) {
							continue;
						}
						BGR& bgr = img.ptr<BGR>(ry + 200)[rx + 200];

						for (int zh = i - 20; zh < i + 20; zh++) {
							if (zh < 0) continue;
							if (zh > 200) break;
							uint16_t green = median594[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * zh];
							if (bgr.green < 255 * 255 && bgr.green < (green - green_pixel_lower) * (255 * 255 / green_pixel_upper)) {
								if (green > green_pixel_upper) bgr.green = 255 * 255;
								else bgr.green = (green - green_pixel_lower) * (255 * 255 / green_pixel_upper);
							}
							uint16_t blue = stack[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * zh];

							if (bgr.blue < 255 * 255 && bgr.blue < (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper)) {
								if (blue > blue_pixel_upper) bgr.blue = 255 * 255;
								else bgr.blue = (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper);
							}

							uint16_t red = median640[(nx + rx) + (ny + ry) * 2048 + 2048 * 2048 * zh];
							if (bgr.red < 255 * 255 && bgr.red < (red - red_pixel_lower) * (255 * 255 / red_pixel_upper)) {
								if (red > red_pixel_upper) bgr.red = 255 * 255;
								else bgr.red = (red - red_pixel_lower) * (255 * 255 / red_pixel_upper);
							}
						}
					}
				}

				// for all the dots, show their outline
				for (int j = 0; j < nuc->close_dots594.size(); j++) {
					int dx, dy, dz;
					std::tie(dx, dy, dz) = nuc->close_dots594.at(j)->local_max;
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(0, 255 * 255, 0), 1, 8, 0);
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);
				}

				// for all the dots, show their outline
				for (int j = 0; j < nuc->close_dots640.size(); j++) {
					int dx, dy, dz;
					std::tie(dx, dy, dz) = nuc->close_dots640.at(j)->local_max;
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(0, 0, 255 * 255), 1, 8, 0);
					cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);
				}

				for (int j = 0; j < nuc->boundary.size(); j++) {
					int x, y, z;

					std::tie(x, y, z) = nuc->boundary.at(j);

					if (x - nx > 200 || x - nx < -200 || y - ny > 200 || y - ny < -200) continue;
					if (z == i) {
						BGR& bgr = img.ptr<BGR>(y - ny + 200)[x - nx + 200];
						bgr.red = 255 * 255;
						bgr.green = 255 * 255;
						bgr.blue = 255 * 255;

					}
				}
				//
				char filename[500];
				sprintf(filename, "%smaxprojects/pos%d_nuc%d.png", dir, pos, k);
				cv::imwrite(filename, img);

			}
			std::cout << "done." << std::endl;

			// make upper and lower limits

			// show up close version of cell, with dots
			// nuclei is dapi
			// 594 is 
			// cell radius is 50
			// so maybe 150 by 150 by 75 box


			// compare to ground truth
			/*
			FILE* groundtruth_405 = fopen("405_groundtruth.csv", "r");
			std::vector<std::tuple<int, int, int>> groundtruth_405_maxima;
			int x, y, z;
			char buf[500];
			fscanf(groundtruth_405, "%s\n", buf);
			while (fscanf(groundtruth_405, "%d,%d,%d\n", &x, &y, &z) != EOF) {
				groundtruth_405_maxima.push_back(std::make_tuple(x, y, z));
			}
			fclose(groundtruth_405);
			std::cout << "done loading groundtruth" << std::endl;

			std::vector<std::tuple<int, int, int>> maxima_notfound;
			int found = 0;
			for (int i = 0; i < nucleii.size(); i++) {
				int max_x, max_y, max_z;
				std::tie(max_x, max_y, max_z) = nucleii.at(i)->local_max;
				int imatch = -1;
				for (int j = 0; j < groundtruth_405_maxima.size(); j++) {
					int true_x, true_y, true_z;
					std::tie(true_x, true_y, true_z) = groundtruth_405_maxima.at(j);
					if ((max_x - true_x) * (max_x - true_x) + (max_y - true_y) * (max_y - true_y) + (max_z - true_z) * (max_z - true_z) < 50 * 50) {
						imatch = j;
					}
				}
				if (imatch == -1) {
					maxima_notfound.push_back(std::make_tuple(max_x, max_y, max_z));
				}
				else {
					found++;
				}
			}
			if (maxima_notfound.size() == 0) {
				std::cout << "All maxima were found in groundtruth." << std::endl;
			}
			else {
				std::cout << "The following " << maxima_notfound.size() << " maxima were not found in the groundtruth : " << std::endl;
				for (int i = 0; i < maxima_notfound.size(); i++) {
					int x, y, z;
					std::tie(x, y, z) = maxima_notfound.at(i);
					std::cout << x << " " << y << " " << z << std::endl;
				}
			}

			std::vector<std::tuple<int, int, int>> groundtruth_405_notfound;
			found = 0;
			for (int i = 0; i < groundtruth_405_maxima.size(); i++) {
				int max_x, max_y, max_z;
				std::tie(max_x, max_y, max_z) = groundtruth_405_maxima.at(i);
				int imatch = -1;
				for (int j = 0; j < nucleii.size(); j++) {
					int true_x, true_y, true_z;
					std::tie(true_x, true_y, true_z) = nucleii.at(j)->local_max;
					if ((max_x - true_x) * (max_x - true_x) + (max_y - true_y) * (max_y - true_y) + (max_z - true_z) * (max_z - true_z) < 50 * 50) {
						imatch = j;
					}
				}
				if (imatch == -1) {
					groundtruth_405_notfound.push_back(std::make_tuple(max_x, max_y, max_z));
				}
				else {
					found++;
				}
			}
			if (groundtruth_405_notfound.size() == 0) {
				std::cout << "All groundtruth were found in maxima." << std::endl;
			}
			else {
				std::cout << "The following " << groundtruth_405_notfound.size() << " groundtruth were not found in the maxima: " << std::endl;
				for (int i = 0; i < groundtruth_405_notfound.size(); i++) {
					int x, y, z;
					std::tie(x, y, z) = groundtruth_405_notfound.at(i);
					std::cout << x << " " << y << " " << z << std::endl;
				}
			}

			// give nuclei their id

			// match dots with their closest nuclei, also keep track of second closest.

			// this is the nuclei, maybe we need an id to go along with it.
			*/
			delete[] stack;
			delete[] filtered;
			delete[] gradientField;
			delete[] stack488;
			delete[] median594;
			delete[] median640;
			for (Nucleus* nuc : nucleii) {
				delete nuc;
			}
			for (Dot* dot : dots594) {
				delete dot;
			}
			for (Dot* dot : dots647) {
				delete dot;
			}
			pos++;
			sprintf(dapi_filename, "%s#%d_New.tif", dir, pos);
		}
	}
};