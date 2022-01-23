#include <iostream>
#include <cmath>
#include <tuple>
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
};


void findMaxima(float* voxels, std::vector<std::tuple<int, int, int>> &maxima) {
	for (int z = 1; z < 200; z++) {
		for (int y = 1; y < 2047; y++) {
			for (int x = 1; x < 2047; x++) {
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
				maxima.push_back(std::make_tuple(x, y, z));
			}
		}
	}
}

void gradientField3d(float* voxels, float3* gradientField) {
	float sixth_order_centered[4] = {0, 3.0 / 4, -3.0 / 20, 1.0 / 60};
	float sixth_order_forward[7] = {-49.0 / 20, 6.0, -15.0 / 2, 20.0 / 3, -15.0 / 4, 6.0 / 5, -1.0 / 6};
	for (int z = 0; z < 201; z++) {
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
						dz += voxels[x + 2048 * y + 2048 * 2048 * (z+k)] * sixth_order_forward[k];
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
	}
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

void segment_blob(
	std::vector<std::tuple<int, int, int>>& points, 
	std::vector<std::tuple<int,int,int>>& boundary, 
	std::tuple<int, int, int>& local_max, 
	float* laplacian,
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
	std::cout << "Computing blob at maximum: (" << x_orig << ", " << y_orig << ", " << z_orig << ")" << std::endl; 
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
			
			float next_laplacian = laplacian[x2 + 2048 * y2 + 2048 * 2048 * z2];

			// instead of using laplacian, we look at whether the rate of change of the radial second derivative is positive.
			// so compute d^2 f / dr^2 at new point, if greater than zero, end.
			// d^2 f / dr^2 = grad ( grad(f) dot r ) dot r
			
			float sixth_order_centered[4] = { 0, 3.0 / 4, -3.0 / 20, 1.0 / 60 };
			float sixth_order_forward[7] = { -49.0 / 20, 6.0, -15.0 / 2, 20.0 / 3, -15.0 / 4, 6.0 / 5, -1.0 / 6 };
			float gradx = 0, grady = 0, gradz = 0;

			float3 r;
			r.x = x2 - x_orig;
			r.y = y2 - y_orig;
			r.z = z2 - z_orig;
			float r_mag = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
			r.x = r.x / r_mag;
			r.y = r.y / r_mag;
			r.z = r.z / r_mag;
			
			if (x2 < 3) {
				for (int j = 0; j < 7; j++) {
					float3 gradf = gradientField[(x2+j) + 2048 * y2 + 2048 * 2048 * z2];
					//float3 r;
					//r.x = x2 + j - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradx += gradfdotr * sixth_order_forward[j];
				}
			}
			else if (x2 > 2048 - 4) {
				for (int j = 0; j< 7; j++) {
					float3 gradf = gradientField[(x2 - j) + 2048 * y2 + 2048 * 2048 * z2];
					//float3 r;
					//r.x = x2 - j - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradx += -gradfdotr * sixth_order_forward[j];
				}
			}
			else {
				for (int j = 1; j < 4; j++) {
					float3 gradf = gradientField[(x2 - j) + 2048 * y2 + 2048 * 2048 * z2];
					//float3 r;
					//r.x = x2 - j - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradx += -gradfdotr * sixth_order_centered[j];

					gradf = gradientField[(x2 + j) + 2048 * y2 + 2048 * 2048 * z2];
					//r.x = x2 + j - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 - z_orig;
					gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradx += gradfdotr * sixth_order_centered[j];
				}
			}


			if (y2 < 3) {
				for (int k = 0; k < 7; k++) {
					float3 gradf = gradientField[x2 + 2048 * (y2+k) + 2048 * 2048 * z2];
					//float3 r;
					//r.x = x2 - x_orig;
					//r.y = y2+k - y_orig;
					//r.z = z2 - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					grady += gradfdotr * sixth_order_forward[k];
				}
			}
			else if (y2 > 2048 - 4) {
				for (int k = 0; k < 7; k++) {
					float3 gradf = gradientField[x2 + 2048 * (y2 - k) + 2048 * 2048 * z2];
					//float3 r;
					//r.x = x2 - x_orig;
					//r.y = y2 - k - y_orig;
					//r.z = z2 - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;// ) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					grady += -gradfdotr * sixth_order_forward[k];
				}
			}
			else {
				for (int k = 1; k < 4; k++) {
					float3 gradf = gradientField[x2 + 2048 * (y2 + k) + 2048 * 2048 * z2];
					//float3 r;
					//r.x = x2 - x_orig;
					//r.y = y2 + k - y_orig;
					//r.z = z2 - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					grady += gradfdotr * sixth_order_centered[k];

					gradf = gradientField[x2 + 2048 * (y2 - k) + 2048 * 2048 * z2];
					//r.x = x2 - x_orig;
					//r.y = y2 - k - y_orig;
					//r.z = z2 - z_orig;
					gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					grady += -gradfdotr * sixth_order_centered[k];
				}
			}


			if (z2 < 3) {

				for (int l = 0; l < 7; l++) {
					float3 gradf = gradientField[x2 + 2048 * y2  + 2048 * 2048 * (z2+l)];
					//float3 r;
					//r.x = x2 - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 + l -z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradz += gradfdotr * sixth_order_forward[l];
				}
			}
			else if (z2 > 201 - 4) {
				for (int l = 0; l < 7; l++) {
					float3 gradf = gradientField[x2 + 2048 * y2 + 2048 * 2048 * (z2 - l)];
					//float3 r;
					//r.x = x2 - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 - l - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradz += -gradfdotr * sixth_order_forward[l];
				}
			}
			else {

				for (int l = 1; l < 4; l++) {
					float3 gradf = gradientField[x2 + 2048 * y2 + 2048 * 2048 * (z2 + l)];

					//float3 r;
					//r.x = x2 - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 + l - z_orig;
					float gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradz += gradfdotr * sixth_order_centered[l];

					gradf = gradientField[x2 + 2048 * y2 + 2048 * 2048 * (z2 - l)];
					//r.x = x2 - x_orig;
					//r.y = y2 - y_orig;
					//r.z = z2 - l - z_orig;
					gradfdotr = gradf.x * r.x + gradf.y * r.y + gradf.z * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					gradz += -gradfdotr * sixth_order_centered[l];

				}
			}
			float3 next_gradient = gradientField[x2 + 2048 * y2 + 2048 * 2048 * z2];
			float d2fdr2 = gradx * r.x + grady * r.y + gradz * r.z;//) / sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
			if (d2fdr2 >= 0 
				//next_laplacian>= 0 // || next_laplacian < this_laplacian //-1/(sigma*sigma)*0.37
				|| next_gradient.x * r.x + next_gradient.y * r.y + next_gradient.z * r.z >= 0
				//|| next_gaussian > this_gaussian // solves kissing problem?
				|| x2 == 0 || y2 == 0 || z2 == 0 || x2 == 2047 || y2 == 2047 || z2 == 200) {
				// if laplacian is greater than or equal 0, then we're in the boundary.
				// also if on boundary of image
				boundary.push_back(std::make_tuple(x2, y2, z2));
			} else {
				// if laplacian is less than 0, then we're in points
				points.push_back(std::make_tuple(x2, y2, z2));
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
float* gaussian_filter3D_parallel(uint16_t* input, int sigmaxy, int sigmaz) {
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

	float* result = new float[2048 * 2048 * 201];
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
	return result;
}


std::vector<cv::Mat*>* load_tiff(const char* filename) {
	TIFF* tiff = TIFFOpen(filename, "r");
	// turn off warnings
	TIFFSetWarningHandler(0);

	std::vector<cv::Mat*>* mats = new std::vector<cv::Mat*>();
	do {
		cv::Mat* loaded = new cv::Mat(2048, 2048, CV_16U);

		for (int i = 0; i < 2048; i++) {
			TIFFReadScanline(tiff, (loaded->data + 4096 * i), i);
		}
		mats->push_back(loaded);
	} while (TIFFReadDirectory(tiff));
	TIFFClose(tiff);
	return mats;
}

// Loading tiff data into a packed array form
// at this point not flexible to dimensions of array.
uint16_t* load_tiff2(const char* filename) {
	TIFF* tiff = TIFFOpen(filename, "r");
	uint16_t* output = new uint16_t[2048 * 2048 * 201];

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
	return output;
}

int main()
{
	int threshold_405_lower = 160;
	int threshold_405_higher = 170;
	//std::cout << "please observe this dog..." << std::endl;
	std::cout << "Loading DAPI...";
	//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
	uint16_t* stack = load_tiff2("endogenous_0x_#1_New.tif");
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
				} else {
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
	int sigma = 25;
	float* gaussian_result = gaussian_filter3D_parallel(filtered, sigma, sigma/2);
	std::cout << "done." << std::endl;

	time(&end);
	std::cout << "gaussian filter took " << difftime(end, start) << " seconds." << std::endl;

	std::vector<std::tuple<int, int, int>> maxima;
	time(&start);
	findMaxima(gaussian_result, maxima);
	time(&end);
	std::cout << "Found " << maxima.size() << " nuclei." << std::endl;
	std::cout << "findMaxima took " << difftime(end, start) << " seconds." << std::endl;

	//std::cout << "Redoing gaussian filter with smaller sigma for edge detection...";
	//gaussian_result = gaussian_filter3D_parallel(filtered, sigma/2, sigma / 4);
	//std::cout << "done." << std::endl;

	float* laplacian = new float[2048 * 2048 * 201];
	std::cout << "computing laplacian...";             
	time(&start);
	laplacianFilter3D(gaussian_result, laplacian);
	std::cout << "done." << std::endl; 
	time(&end);
	std::cout << "Laplacian filter took: " << difftime(end, start) << " seconds." << std::endl;

	int gradientfieldsize = 2048 * 2048 * 201;
	float3* gradientField = new float3[gradientfieldsize];
	std::cout << "computing gradient field...";
	time(&start);
	gradientField3d(gaussian_result, gradientField);
	std::cout << "done." << std::endl;
	time(&end);
	std::cout << "Gradient field took: " << difftime(end, start) << " seconds." << std::endl;

	//std::cout << "Testing gradient field" << std::endl;
	//std::cout << "test test: " << gradientField[2047 + 2048 * 2047 + 2048 * 2048 * 200].y << std::endl;
	//return 0;

	std::cout << "segmenting blobs...";
	std::vector<Blob*> blobs;
	time(&start);
	for (int j = 0; j < maxima.size(); j++) {
		Blob* blob = new Blob;
		blob->points.reserve(200000);
		blob->boundary.reserve(20000);
		blob->local_max = maxima.at(j);
		segment_blob(blob->points, blob->boundary, maxima.at(j), laplacian, gradientField);
		blobs.push_back(blob);
		std::cout << "blob has " << blob->points.size() << " voxels" << std::endl;
		std::cout << "boundary has " << blob->boundary.size() << " voxels" << std::endl;
	}
	time(&end);

	std::cout << "done." << std::endl;
	std::cout << "blob segmentation took: " << difftime(end, start) << " seconds." << std::endl;

	//std::cout << "please observe this dog..." << std::endl;
	std::cout << "Loading 594...";
	//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
	uint16_t* stack594 = load_tiff2("endogenous_0x_#1_New_2.tif");
	std::cout << " done." << std::endl;

	pointer = stack594;
	for (int z = 0; z < 201; z++) {
		for (int y = 0; y < 2048; y++) {
			for (int x = 0; x < 2048; x++) {
				if ((*pointer) <= 190) {
					*pointer = 0;
				}
				else {
					*pointer = (*pointer) - 190;
				}
				pointer++;
			}
		}
	}
	std::cout << "computing 594 gaussian...";
	float* gaussian_594 = gaussian_filter3D_parallel(stack594, 5, 2);
	std::cout << "done." << std::endl; 

	std::vector<std::tuple<int, int, int>> maxima594;
	findMaxima(gaussian_594, maxima594);
	std::cout << "Found " << maxima594.size() << " dots." << std::endl;

	//std::cout << "please observe this dog..." << std::endl;
	std::cout << "Loading 647...";
	//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
	uint16_t* stack647 = load_tiff2("endogenous_0x_#1_New_4.tif");
	std::cout << " done." << std::endl;
	pointer = stack647;
	for (int z = 0; z < 201; z++) {
		for (int y = 0; y < 2048; y++) {
			for (int x = 0; x < 2048; x++) {
				if ((*pointer) <= 190) {
					*pointer = 0;
				}
				else {
					*pointer = (*pointer) - 190;
				}
				pointer++;
			}
		}
	}
	std::cout << "computing 647 gaussian...";
	float* gaussian_647 = gaussian_filter3D_parallel(stack647, 5, 2);
	std::cout << "done." << std::endl;

	std::vector<std::tuple<int, int, int>> maxima647;
	findMaxima(gaussian_647, maxima647);
	std::cout << "Found " << maxima647.size() << " dots." << std::endl;

	int i = 0;
	int k = 0;
	bool blink = true;
	int key = 0;
	do {
		//cv::Mat img(2048, 2048, CV_32F, gaussian_result + 2048*2048*i);
		cv::Mat img(2048, 2048, CV_16U, stack + 2048 * 2048 * i);
		cv::Mat dst(2048, 2048, CV_16UC3);
		cv::cvtColor(img, dst, cv::COLOR_GRAY2BGR);

		for (int j = 0; j < maxima.size(); j++) {
			int z = std::get<2>(maxima.at(j));
			if (abs(i - z) < 10) {
				int cx = std::get<0>(maxima.at(j));
				int cy = std::get<1>(maxima.at(j));
				cv::rectangle(dst,cv::Point(cx - 5, cy - 5), cv::Point(cx + 5, cy + 5), cv::Scalar(0, 0, 65535), cv::FILLED, cv::LINE_8);
			}
		}
		if (blink && blobs.size() > 0) {
			Blob* blob = blobs.at(k);

			for (int j = 0; j < blob->boundary.size(); j++) {
				int x, y, z;

				std::tie(x, y, z) = blob->boundary.at(j);

				if (z == i) {
					BGR& bgr = dst.ptr<BGR>(y)[x];
					bgr.red = 0;
					bgr.green = 255;
					bgr.blue = 0;
				}
			}
		}
		//cv::Mat resized(512, 512, CV_32F);
		cv::Mat resized(512, 512, CV_16UC3);
		cv::resize(dst, resized, cv::Size(512, 512));
		//cv::imshow("Display window", resized);
		cv::imshow("Display window", resized*255);
		key = cv::waitKey(0);
		if (key == 's') {
			i = (i + 1) % 201;
		}
		else if (key == 'w') {
			i = (201 + i - 1) % 201;
		}
		else if (key == 'a') {
			k = (k + 1) % blobs.size();
			i = std::get<2>(blobs.at(k)->local_max);
		}
		else if (key == 'd') {
			k = (blobs.size() + k - 1) % blobs.size();
			i = std::get<2>(blobs.at(k)->local_max);
		}
		else {
			//blink = !blink;
		}

	} while (key != 27);

	i = 0;
	k = 0;
	do {
		cv::Mat img(2048, 2048, CV_32F, laplacian + 2048 * 2048 * i);
		double min, max;
		cv::minMaxLoc(img, &min, &max);
		img = (img - min) / (max - min);

		cv::Mat dst(2048, 2048, CV_32FC3);
		cv::cvtColor(img, dst, cv::COLOR_GRAY2BGR);

		//cv::Mat img(2048, 2048, CV_16U, stack + 2048 * 2048 * i);

		for (int j = 0; j < maxima.size(); j++) {
			int z = std::get<2>(maxima.at(j));
			if (abs(i - z) < 10) {
				int cx = std::get<0>(maxima.at(j));
				int cy = std::get<1>(maxima.at(j));
				cv::rectangle(dst, cv::Point(cx - 5, cy - 5), cv::Point(cx + 5, cy + 5), cv::Scalar(0, 0, 1), cv::FILLED, cv::LINE_8);
			}
		}
		Blob* blob = blobs.at(k);

		for (int j = 0; j < blob->boundary.size(); j++) {
			int x, y, z;

			std::tie(x, y, z) = blob->boundary.at(j);

			if (z == i) {
				BGR_float& bgr = dst.ptr<BGR_float>(y)[x];
				bgr.red = 0;
				bgr.green = 1;
				bgr.blue = 0;
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
			k = (k + 1) % blobs.size();
			i = std::get<2>(blobs.at(k)->local_max);
		}
		else if (key == 'd') {
			k = (blobs.size() + k - 1) % blobs.size();
			i = std::get<2>(blobs.at(k)->local_max);
		}
		else {
			//blink = !blink;
		}

	} while (key != 27);

	// compare to ground truth
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
	for (int i = 0; i < maxima.size(); i++) {
		int max_x, max_y, max_z;
		std::tie(max_x, max_y, max_z) = maxima.at(i);
		int imatch = -1;
		for (int j = 0; j < groundtruth_405_maxima.size(); j++) {
			int true_x, true_y, true_z; 
			std::tie(true_x, true_y, true_z) = groundtruth_405_maxima.at(j);
			if ((max_x - true_x) * (max_x - true_x) + (max_y - true_y) * (max_y - true_y) + (max_z - true_z) * (max_z - true_z)  < 50 * 50) {
				imatch = j;
			}
		}
		if (imatch == -1) {
			maxima_notfound.push_back(std::make_tuple(max_x, max_y, max_z));
		} else {
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
		for (int j = 0; j < maxima.size(); j++) {
			int true_x, true_y, true_z;
			std::tie(true_x, true_y, true_z) = maxima.at(j);
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

	FILE* maxima_csv = fopen("405_maxima.csv", "w");
	fprintf(maxima_csv, "x,y,z\n");
	for (int i = 0; i < maxima.size(); i++) {
		int max_x, max_y, max_z;
		std::tie(max_x, max_y, max_z) = maxima.at(i);
		fprintf(maxima_csv, "%d,%d,%d\n", max_x, max_y, max_z);
	}
	fclose(maxima_csv);
};