#pragma once

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


#define GAUSSIAN_CUTOFF 6
#define M_PI 3.14159265358979323846


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
	const int size_upper_limit;
	const int size_lower_limit;
public:
	Nucleus(int sul, int sll) : size_upper_limit(sul), size_lower_limit(sll) {
	}

	int id;
	std::vector<Dot*> close_dots594;
	std::vector<Dot*> close_dots640;
	std::vector<std::tuple<int, int, int>> close_points594;
	std::vector<std::tuple<int, int, int>> close_points640;


	bool validSize() {
		return points.size() >= size_lower_limit && points.size() <= size_upper_limit;
	}
};

class Dot : public Blob {
	const int size_upper_limit;
	const int size_lower_limit;
public:
	int id;

	Dot(const int sul, const int sll) : size_upper_limit(sul), size_lower_limit(sll) {
	}

	bool validSize() {
		return points.size() >= size_lower_limit && points.size() <= size_upper_limit;
	}
};



int deltaXSq(std::tuple<int, int, int> X, std::tuple<int, int, int> Y) {
	int x1, y1, z1, x2, y2, z2;
	std::tie(x1, y1, z1) = X;
	std::tie(x2, y2, z2) = Y;
	x1 = x1 - x2;
	y1 = y1 - y2;
	z1 = z1 - z2;
	return x1 * x1 + y1 * y1 + z1 * z1;
}
float deltaX(std::tuple<int, int, int> X, std::tuple<int, int, int> Y) {
	return sqrt(deltaXSq(X, Y));

}


void getTiffDimensions(const char* filename, int& width, int& height, int& depth) {
	TIFF* tiff = TIFFOpen(filename, "r");
	TIFFSetWarningHandler(0);


	TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);

	depth = 0;
	while (TIFFSetDirectory(tiff, depth)) {
		depth++;
	}

	TIFFClose(tiff);
}

void loadTiff(uint16_t* output, const char* filename, const int width, const int height, const int depth) {
		TIFF* tiff = TIFFOpen(filename, "r");
		// turn off warnings
		TIFFSetWarningHandler(0);
		int zed = 0, z = 0;
		do {
			int ied = 0;
			for (int i = 0; i < height; i++) {
				TIFFReadScanline(tiff, (output + ied + zed), i);
				ied += width;
			}
			zed += width * height;
			z++;
		} while (TIFFReadDirectory(tiff) && z < depth);
		TIFFClose(tiff);
}


void findMaxima(float* voxels, const int width, const int height, const int depth, std::vector<std::tuple<int, int, int>>& maxima) {
	for (int z = 2; z < depth-2; z++) {
		for (int y = 2; y < height -2; y++) {
			for (int x = 2; x < width - 2; x++) {
				float current = voxels[width * y + height* width * z + x];
				// 6 checks
				if (current <= voxels[width * y + height* width * z + x + 1])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1])
					continue;
				if (current <= voxels[width * y + height* width * z + x + width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 1 + width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 - width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 - width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 1 + width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 - width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 - width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + width + width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + width - width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - width + width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - width - width*height])
					continue;

				// go out to distance of 2 to avoid numerical error

				// x - 2 face
				if (current <= voxels[width * y + height* width * z + x - 2 - 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 2 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 2 - 1 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 1 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 1 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 1 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 1 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 2 - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 2 + 1 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 1 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 1 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 1 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 1 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 2 + 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 + 2 * width + 2 * width*height])
					continue;

				// x + 2 face
				if (current <= voxels[width * y + height* width * z + x + 2 - 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 2 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 2 - 1 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 1 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 1 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 1 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 1 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 2 - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 2 + 1 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 1 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 1 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 1 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 1 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 2 + 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 + 2 * width + 2 * width*height])
					continue;


				// y - 2 face
				if (current <= voxels[width * y + height* width * z + x - 1 - 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 - 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 - 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 - 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 - 2 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 1 - 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 - 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 - 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 - 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 - 2 * width + 2 * width*height])
					continue;


				// y + 2 face
				if (current <= voxels[width * y + height* width * z + x - 1 + 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + 2 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 1 + 2 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + 2 * width - 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + 2 * width])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + 2 * width + 1 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + 2 * width + 2 * width*height])
					continue;

				// z - 2 face
				if (current <= voxels[width * y + height* width * z + x - 1 - 1 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + 1 * width - 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 1 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 * width - 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 1 - 1 * width - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 - 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + 1 * width - 2 * width*height])
					continue;


				// z + 2 face
				if (current <= voxels[width * y + height* width * z + x - 1 - 1 * width + 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x - 1 + 1 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x - 1 * width + 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 * width + 2 * width*height])
					continue;

				if (current <= voxels[width * y + height* width * z + x + 1 - 1 * width + 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + 2 * width*height])
					continue;
				if (current <= voxels[width * y + height* width * z + x + 1 + 1 * width + 2 * width*height])
					continue;

				maxima.push_back(std::make_tuple(x, y, z));
			}
		}
	}
}


void gradientField3d(float* voxels, const int width, const int height, const int depth, float3* gradientField) {
	float sixth_order_centered[4] = { 0, 3.0 / 4, -3.0 / 20, 1.0 / 60 };
	float sixth_order_forward[7] = { -49.0 / 20, 6.0, -15.0 / 2, 20.0 / 3, -15.0 / 4, 6.0 / 5, -1.0 / 6 };
	concurrency::parallel_for(0, depth, [&sixth_order_centered, &sixth_order_forward, gradientField, voxels, width, height, depth](int z) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				float dx = 0, dy = 0, dz = 0;

				if (x < 3) {
					for (int i = 0; i < 7; i++) {
						dx += voxels[(x + i) + width * y + width * height * z] * sixth_order_forward[i];
					}
				}
				else if (x > width - 4) {
					for (int i = 0; i < 7; i++) {
						dx += -voxels[(x - i) + width * y + width * height * z] * sixth_order_forward[i];
					}
				}
				else {
					for (int i = 1; i < 4; i++) {
						dx += voxels[(x + i) + width * y + width * height * z] * sixth_order_centered[i];
						dx += -voxels[(x - i) + width * y + width * height * z] * sixth_order_centered[i];
					}
				}
				gradientField[x + width * y + width * height * z].x = dx;

				if (y < 3) {
					for (int j = 0; j < 7; j++) {
						dy += voxels[x + width * (y + j) + width * height * z] * sixth_order_forward[j];
					}
				}
				else if (y > height - 4) {
					for (int j = 0; j < 7; j++) {
						dy += -voxels[x + width * (y - j) + width * height * z] * sixth_order_forward[j];
					}
				}
				else {
					for (int j = 1; j < 4; j++) {
						dy += voxels[x + width * (y + j) + width * height * z] * sixth_order_centered[j];
						dy += -voxels[x + width * (y - j) + width * height * z] * sixth_order_centered[j];
					}
				}
				gradientField[x + width * y + width * height * z].y = dy;

				if (z < 3) {
					for (int k = 0; k < 7; k++) {
						dz += voxels[x + width * y + width * height * (z + k)] * sixth_order_forward[k];
					}
				}
				else if (z > depth - 4) {
					for (int k = 0; k < 7; k++) {
						dz += -voxels[x + width * y + width * height * (z - k)] * sixth_order_forward[k];
					}
				}
				else {
					for (int k = 1; k < 4; k++) {
						dz += voxels[x + width * y + width * height * (z + k)] * sixth_order_centered[k];
						dz += -voxels[x + width * y + width * height * (z - k)] * sixth_order_centered[k];

					}
				}
				gradientField[x + width * y + width * height * z].z = dz;
			}
		}
		});
}


void eigenvalues(std::array<std::array<float, 3>, 3> hessian, float& eig1, float& eig2, float& eig3) {
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
		if (hessian[i][i] == 0) {
			for (int j = i + 1; j < 3; j++) {
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
				if (j != i) tmp[i] -= hessian[i][j] / hessian[i][i];
			}
		}
	}

	float ev_mag = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);

	ev.x = tmp[0] / ev_mag;
	ev.y = tmp[1] / ev_mag;
	ev.z = tmp[2] / ev_mag;
}


// medium fast gaussian filter
// In this function, for optimization purposes I try to use some pointer
// arithmetic so that I'm always reading from continuous memory.
template <typename T>
void gaussian_filter3D_parallel(T* input, const int width, const int height, const int depth, int sigmaxy, int sigmaz, float* result) {
	float float_sigmaxy = (float)sigmaxy;
	float float_sigmaz = (float)sigmaz;
	// prerun expensive exp operation in 1d array
	float* kernelxy = new float[2 * GAUSSIAN_CUTOFF * sigmaxy];
	float* kernelz = new float[2 * GAUSSIAN_CUTOFF * sigmaz];
	float normxy = 1.0 / sqrt(2 * M_PI * float_sigmaxy * float_sigmaxy);
	float normz = 1.0 / sqrt(2 * M_PI * float_sigmaz * float_sigmaz);

	for (int r = 0; r <= 2 * GAUSSIAN_CUTOFF * sigmaxy; r++) {
		kernelxy[r] = exp(-r * r / (2 * float_sigmaxy * float_sigmaxy)) * normxy;
	}
	for (int r = 0; r <= 2 * GAUSSIAN_CUTOFF * sigmaz; r++) {
		kernelz[r] = exp(-r * r / (2 * float_sigmaz * float_sigmaz)) * normz;
	}

	float* temp = new float[width * height * depth];

	concurrency::static_partitioner partitioner;
	//concurrency::critical_section cs;
	uint16_t* yaxis = new uint16_t[height*width];
	float* zaxis = new float[depth*height];
	float* xaxis = new float[width*depth]; 

	concurrency::parallel_for(0, width, [&input, &result, &sigmaxy, &kernelxy, &temp, &yaxis, depth, width, height](int x) {
		for (int z = 0; z < depth; z++) {
			// copy to contiguous memory for better speed.
			//const  int h = height;

			for (int y = 0; y < height; y++) {
				yaxis[height * x + y] = input[x + width * height * z + y * width];
				//result[z * height + height * depth * x + y] = input[x + width * height * z + y * width];
			}
			//continue;

			for (int y = 0; y < height; y++) {
				float sum = 0;
				int j_min, j_max;
				if (y - GAUSSIAN_CUTOFF * sigmaxy > 0) {
					j_min = y - GAUSSIAN_CUTOFF * sigmaxy;
				}
				else {
					j_min = 0;
				}
				if (y + GAUSSIAN_CUTOFF * sigmaxy < height) {
					j_max = y + GAUSSIAN_CUTOFF * sigmaxy;
				}
				else {
					j_max = height;
				}

				for (int j = j_min; j < y; j++) {
					sum += yaxis[height*x + j] * kernelxy[y - j];
				}
				for (int j = y; j < j_max; j++) {
					sum += yaxis[height*x + j] * kernelxy[j - y];
				}
				result[z * height + height * depth * x + y] = sum;
			}
		}
		}, partitioner);
	//}

	concurrency::parallel_for(0, height, [&result, &temp, &sigmaz, &kernelz, &zaxis, width, height, depth](int y) {
		for (int x = 0; x < width; x++) {
			// copy to contiguous memory
			for (int z = 0; z < depth;z++) {
				zaxis[depth*y + z] = result[y + height * depth * x + z * height];
				//temp[depth * x + width * depth * y + z] = result[y + height * depth * x + z * height];
			}
			//continue;
			for (int z = 0; z < depth; z++) {
				int j_min, j_max;
				float sum = 0;
				if (z - GAUSSIAN_CUTOFF * sigmaz > 0) {
					j_min = z - GAUSSIAN_CUTOFF * sigmaz;
				}
				else {
					j_min = 0;
				}
				if (z + GAUSSIAN_CUTOFF * sigmaz < depth) {
					j_max = z + GAUSSIAN_CUTOFF * sigmaz;
				}
				else {
					j_max = depth;
				}

				for (int j = j_min; j < z; j++) {
					sum += zaxis[depth*y + j] * kernelz[z - j];
				}
				for (int j = z; j < j_max; j++) {
					sum += zaxis[depth*y + j] * kernelz[j - z];
				}
				temp[depth * x + width * depth * y + z] = sum;
			}
		}
		}, partitioner);
	//}

	concurrency::parallel_for(0, depth, [&temp, &result, &sigmaxy, &kernelxy, &xaxis, width, height, depth](int z) {
		for (int y = 0; y < height; y++) {
			int j_min, j_max;
			for (int x = 0; x < width; x++) {
				xaxis[width*z + x] = temp[z + width * depth * y + x * depth];
				//result[width * y + width * height * z + x] = temp[z + width * depth * y + x * depth];
			}
			//continue;
			for (int x = 0; x < width; x++) {
				float sum = 0;
				if (x - GAUSSIAN_CUTOFF * sigmaxy > 0) {
					j_min = x - GAUSSIAN_CUTOFF * sigmaxy;
				}
				else {
					j_min = 0;
				}
				if (x + GAUSSIAN_CUTOFF * sigmaxy < width) {
					j_max = x + GAUSSIAN_CUTOFF * sigmaxy;
				}
				else {
					j_max = width;
				}
				for (int j = j_min; j < x; j++) {
					sum += xaxis[width*z + j] * kernelxy[x - j];
				}
				for (int j = x; j < j_max; j++) {
					sum += xaxis[width*z + j] * kernelxy[j - x];
				}
				//cs.lock();
				result[width * y + width * height * z + x] = sum;
				//cs.unlock();
			}
		}
		}, partitioner);

	// clean up
	// delete 
	delete[] zaxis;
	delete[] xaxis;
	delete[] yaxis;
	delete[] temp;
	delete[] kernelxy;
	delete[] kernelz;
}

uint16_t median(uint16_t* arr, int len) {
	uint16_t* A = new uint16_t[len];
	uint16_t* B = new uint16_t[len];

	uint16_t* current;
	uint16_t* next;
	current = A;
	next = B;
	for (int i = 0; i < len; i++) {
		A[i] = arr[i];
	}

	// find nth member of the list
	int start_index = 0;
	int n = len/2;
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
			return pivot;
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

	delete[] A;
	delete[] B;
}



void medianFilter3x3(uint16_t* voxels, const int width, const int height, const int depth, uint16_t* filtered) {
	concurrency::parallel_for(0, depth, [&voxels, &filtered, width, height, depth](int z) {
		for (int y = 0; y < height; y++) {
			for (int x = 1; x < width; x++) {
				if (z == 0 || z == depth - 1 || y == 0 || y == height - 1 || x == 0 || x == width - 1) {
					filtered[x + width * y + width * height * z] = 0;
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
							current[(i + 1) + 3 * (j + 1) + 3 * 3 * (k + 1)] = voxels[(x + k) + width * (y + j) + width * height * (z + i)];
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
						filtered[x + width * y + width * height * z] = pivot;
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


// does not compute laplacian on the boundary
// can be implemented later with a one-sided stencil. 
// but I could not find the terms using a quick google search. 
void laplacianFilter3D(float* voxels, const int width, const int height, const int depth, float* laplacian) {
	// z x y laplacian
	float sixth_order_centered[4] = { -49.0 / 18, 3.0 / 2, -3.0 / 20, 1.0 / 90 };
	float sixth_order_forward[8] = { 469.0 / 90, -223.0 / 10, 879.0 / 20, -949.0 / 18, 41.0, -201.0 / 10, 1019.0 / 180, -7.0 / 10 };

	for (int z = 0; z < depth; z++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				float lp = 0;
				// x
				if (x < 3) {
					for (int i = 0; i < 8; i++) {
						lp += voxels[(x + i) + width * y + width * height * z] * sixth_order_forward[i];
					}
				}
				else if (x > width - 4) {
					for (int i = 0; i < 8; i++) {
						lp += voxels[(x - i) + width * y + width * height * z] * sixth_order_forward[i];
					}
				}
				else {
					lp += voxels[x + width * y + width * height * z] * sixth_order_centered[0];
					for (int i = 1; i < 4; i++) {
						lp += voxels[(x + i) + width * y + width * height * z] * sixth_order_centered[i];
						lp += voxels[(x - i) + width * y + width * height * z] * sixth_order_centered[i];
					}
				}

				// y
				if (y < 3) {
					for (int j = 0; j < 8; j++) {
						lp += voxels[x + width * (y + j) + width * height * z] * sixth_order_forward[j];
					}
				}
				else if (y > height - 4) {
					for (int j = 0; j < 8; j++) {
						lp += voxels[x + width * (y - j) + width * height * z] * sixth_order_forward[j];
					}
				}
				else {
					lp += voxels[x + width * y + width * height * z] * sixth_order_centered[0];
					for (int j = 1; j < 4; j++) {
						lp += voxels[x + width * (y + j) + width * height * z] * sixth_order_centered[j];
						lp += voxels[x + width * (y - j) + width * height * z] * sixth_order_centered[j];

					}
				}

				// z 
				if (z < 3) {
					for (int k = 0; k < 8; k++) {
						lp += voxels[x + width * y + width * height * (z + k)] * sixth_order_forward[k];

					}
				}
				else if (z > depth - 4) {
					for (int k = 0; k < 8; k++) {
						lp += voxels[x + width * y + width * height * (z - k)] * sixth_order_forward[k];

					}
				}
				else {
					lp += voxels[x + width * y + width * height * z] * sixth_order_centered[0];
					for (int k = 1; k < 4; k++) {
						lp += voxels[x + width * y + width * height * (z + k)] * sixth_order_centered[k];
						lp += voxels[x + width * y + width * height * (z - k)] * sixth_order_centered[k];
					}
				}
				laplacian[x + width * y + width * height * z] = lp;
			}
		}
	}
}

void hessianAt(float3* gradientField, const int width, const int height, const int depth,  int x2, int y2, int z2, std::array<std::array<float, 3>, 3>& hessian) {
	float sixth_order_centered[4] = { 0, 3.0 / 4, -3.0 / 20, 1.0 / 60 };
	float sixth_order_forward[7] = { -49.0 / 20, 6.0, -15.0 / 2, 20.0 / 3, -15.0 / 4, 6.0 / 5, -1.0 / 6 };

	if (x2 < 3) {
		for (int j = 0; j < 7; j++) {
			float3 gradf = gradientField[(x2 + j) + width * y2 + width * height * z2];
			hessian[0][0] += gradf.x * sixth_order_forward[j];// d2fdx2
			hessian[1][0] += gradf.y * sixth_order_forward[j]; // d2fdxdy
			hessian[2][0] += gradf.z * sixth_order_forward[j]; // d2fdxdz
		}
	}
	else if (x2 > width - 4) {
		for (int j = 0; j < 7; j++) {
			float3 gradf = gradientField[(x2 - j) + width * y2 + width * height * z2];
			hessian[0][0] -= gradf.x * sixth_order_forward[j];// d2fdx2
			hessian[1][0] -= gradf.y * sixth_order_forward[j]; // d2fdxdy
			hessian[2][0] -= gradf.z * sixth_order_forward[j]; // d2fdxdz
		}
	}
	else {
		for (int j = 1; j < 4; j++) {
			float3 gradf = gradientField[(x2 - j) + width * y2 + width * height * z2];
			hessian[0][0] -= gradf.x * sixth_order_centered[j];// d2fdx2
			hessian[1][0] -= gradf.y * sixth_order_centered[j]; // d2fdxdy
			hessian[2][0] -= gradf.z * sixth_order_centered[j]; // d2fdxdz


			gradf = gradientField[(x2 + j) + width * y2 + width * height * z2];
			hessian[0][0] += gradf.x * sixth_order_centered[j];// d2fdx2
			hessian[1][0] += gradf.y * sixth_order_centered[j]; // d2fdxdy
			hessian[2][0] += gradf.z * sixth_order_centered[j]; // d2fdxdz
		}
	}


	if (y2 < 3) {
		for (int k = 0; k < 7; k++) {
			float3 gradf = gradientField[x2 + width * (y2 + k) + width * height * z2];
			hessian[1][1] += gradf.y * sixth_order_forward[k];// d2fdy2
			hessian[2][1] += gradf.z * sixth_order_forward[k];// d2fdydz
		}
	}
	else if (y2 > height - 4) {
		for (int k = 0; k < 7; k++) {
			float3 gradf = gradientField[x2 + width * (y2 - k) + width * height * z2];
			hessian[1][1] -= gradf.y * sixth_order_forward[k];// d2fdy2
			hessian[2][1] -= gradf.z * sixth_order_forward[k];// d2fdy2
		}
	}
	else {
		for (int k = 1; k < 4; k++) {
			float3 gradf = gradientField[x2 + width * (y2 + k) + width * height * z2];
			hessian[1][1] += gradf.y * sixth_order_centered[k];// d2fdy2
			hessian[2][1] += gradf.z * sixth_order_centered[k];// d2fdy2


			gradf = gradientField[x2 + width * (y2 - k) + width * height * z2];
			hessian[1][1] -= gradf.y * sixth_order_centered[k];// d2fdy2
			hessian[2][1] -= gradf.z * sixth_order_centered[k];// d2fdy2
		}
	}


	if (z2 < 3) {

		for (int l = 0; l < 7; l++) {
			float3 gradf = gradientField[x2 + width * y2 + width * height * (z2 + l)];
			hessian[2][2] += gradf.z * sixth_order_forward[l];// d2fdz2

		}
	}
	else if (z2 > depth - 4) {
		for (int l = 0; l < 7; l++) {
			float3 gradf = gradientField[x2 + width * y2 + width * height * (z2 - l)];
			hessian[2][2] -= gradf.z * sixth_order_forward[l];// d2fdz2
		}
	}
	else {

		for (int l = 1; l < 4; l++) {
			float3 gradf = gradientField[x2 + width * y2 + width * height * (z2 + l)];
			hessian[2][2] += gradf.z * sixth_order_centered[l];// d2fdz2

			gradf = gradientField[x2 + width * y2 + width * height * (z2 - l)];
			hessian[2][2] -= gradf.z * sixth_order_centered[l];// d2fdz2

		}
	}
	hessian[0][1] = hessian[1][0];
	hessian[0][2] = hessian[2][0];
	hessian[1][2] = hessian[2][1];
}



void segment_blob(
	std::vector<std::tuple<int, int, int>>& points,
	std::vector<std::tuple<int, int, int>>& boundary,
	std::tuple<int, int, int>& local_max,
	float3* gradientField,
	const int width,
	const int height,
	const int depth) {
	int back = 0;
	if (!points.size() == 0) {
		std::cout << "starting points size is not 0. Returning out of segment_nucleus" << std::endl;
		return;
	}
	char* visited = new char[width * height * depth];
	memset(visited, 0, sizeof(char) * width * height * depth);
	int x_orig, y_orig, z_orig;
	std::tie(x_orig, y_orig, z_orig) = local_max;
	//std::cout << "Computing blob at maximum: (" << x_orig << ", " << y_orig << ", " << z_orig << ")" << std::endl; 
	//if (max_laplacian > 0) { 
		//std::cout << "MAX LAPLACIAN POSITIVE: " << max_laplacian << std::endl;
	//}
	//else {
		//std::cout << "Max Laplacian negative: " << max_laplacian << std::endl;
	//}
	// breadth first search from the local max
	// radial derivative  d^2 f / dr^2 
	points.push_back(std::make_tuple(x_orig, y_orig, z_orig));

	while (points.size() > back) {
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

			//if (x2 - x_orig > 1023 || y2 - y_orig > 1023 || z2 - z_orig > 100 || x2 - x_orig < -1023 || y2 - y_orig < -1023 || z2 - z_orig < -100) {
				//std::cout << "Error in segment_blob: cell too big: " << x2 - x_orig << " " << y2 - y_orig << " " << z2 - z_orig << std::endl;
				//delete[] visited;
				//return;
			//}

			if (visited[x2 + y2 * width + z2 * width * height] == 1) {
				continue;
			}
			if (x2 == 0 || y2 == 0 || z2 == 0 || x2 == width - 1  || y2 == height - 1 || z2 == depth - 1 ) {
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
				hessianAt(gradientField, width, height, depth, x2, y2, z2, hessian);

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


				float3 next_gradient = gradientField[x2 + height * y2 + height * width * z2];
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
					(eig1 >= 0) || (eig2 >= 0) || (eig3 >= 0)
					) {
					//|| d2fdr2 >= 0 || lambda1 >= 0 || lambda2 >= 0
					//d2fdg2 >=0
					//((eig1 >=0) &&  (eig2>= 0)) || ((eig1>=0) && (eig3>=0)) || ((eig2>=0)&&(eig3>=0)) 
					//eig1+eig2+eig3 >= 0
					//|| d2fdr2 >= 0 
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
			visited[x2 + y2 * width + z2 * width * height] = 1;
		}
	}
	delete[] visited;
}

void scanfloat(float* image, const int width, const int height, const int depth, const float min, const float max) {
	int i = 0;
	int k = 0;
	int key = 0;
	do {
		cv::Mat img(height, width, CV_32F, image + width * height * i);
		cv::Mat resized(512, 512, CV_32F);

		cv::resize(img, resized, cv::Size(512, 512));
		resized = (resized - min) / (max - min);
		cv::imshow("Display window", resized);
		key = cv::waitKey(0);
		if (key == 's') {
			i = (i + 1) % depth;
		}
		else if (key == 'w') {
			i = (depth + i - 1) % depth;
		} 
	} while (key != 27);
}


void scanInt(uint16_t* image, const int width, const int height, const int depth) {
	int i = 0;
	int k = 0;
	int key = 0;
	do {
		cv::Mat img(height, width, CV_16U, image + width * height * i);
		cv::Mat resized(512, 512, CV_16U);
		cv::resize(img, resized, cv::Size(512, 512));
		resized = resized * 255;

		cv::imshow("Display window", resized);
		key = cv::waitKey(0);
		if (key == 's') {
			i = (i + 1) % depth;
		}
		else if (key == 'w') {
			i = (depth + i - 1) % depth;
		} 
	} while (key != 27);
}


void displayBlobsInt(uint16_t* stack, int width, int height, int depth, std::vector<Nucleus*>& nucleii) {
	int i = 0;
	int k = 0;
	bool blink = true;
	int key = 0;
	do {
		cv::Mat img(height, width, CV_16U, stack + width * height * i);
		cv::Mat dst(height, width, CV_16UC3);
		cv::cvtColor(img, dst, cv::COLOR_GRAY2BGR);

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
			i = (i + 1) % depth;
		}
		else if (key == 'w') {
			i = (depth + i - 1) % depth;
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

	} while (key != 27);
}


void displayBlobsFloat(float* stack, int width, int height, int depth, std::vector<Nucleus*>& nucleii) {
	int i = 0;
	int k = 0;
	bool blink = true;
	int key = 0;
	do {
		cv::Mat img(height, width, CV_32F, stack + width * height * i);
		cv::Mat dst(height, width, CV_32FC3);
		cv::cvtColor(img, dst, cv::COLOR_GRAY2BGR);

		for (int l = 0; l < nucleii.size(); l++) {
			Nucleus* nuc = nucleii.at(l);
			std::tuple<int, int, int> max = nuc->local_max;

			if (nuc->validSize()) {
				int z = std::get<2>(max);
				if (abs(i - z) < 10) {
					int cx = std::get<0>(max);
					int cy = std::get<1>(max);
					cv::rectangle(dst, cv::Point(cx - 5, cy - 5), cv::Point(cx + 5, cy + 5), cv::Scalar(0, 0, 1), cv::FILLED, cv::LINE_8);
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
		//cv::Mat resized(512, 512, CV_32F);
		cv::Mat resized(512, 512, CV_32FC3);
		cv::resize(dst, resized, cv::Size(512, 512));
		//cv::imshow("Display window", resized);
		cv::imshow("Display window", resized);
		key = cv::waitKey(0);
		if (key == 's') {
			i = (i + 1) % depth;
		}
		else if (key == 'w') {
			i = (depth + i - 1) % depth;
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

	} while (key != 27);
}

void segmentNuclei(uint16_t* stack, 
	uint16_t* filtered,
	float* gaussian_result,
	float3* gradientField,
	std::vector<Nucleus*>& nuclei, 
	std::vector<std::tuple<int, int, int>>& maxima,
	int threshold_405_lower,
	int threshold_405_higher,
	int sigmaxy, 
	int sigmaz,
	int width, int height, int depth) {

	time_t start, end;
	std::cout << "Computing median filter...";
	time(&start);
	medianFilter3x3(stack, width, height, depth, filtered);
	time(&end);
	std::cout << "done." << std::endl;
	std::cout << "median filter took " << end - start << "seconds" << std::endl;

	uint16_t* pointer = filtered;
	for (int z = 0; z < depth; z++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
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

	std::cout << "Computing gaussian filter...";
	time(&start);
	gaussian_filter3D_parallel<uint16_t>(filtered, width, height, depth, sigmaxy, sigmaz, gaussian_result);
	std::cout << "done." << std::endl;
	time(&end);

	time(&start);
	findMaxima(gaussian_result, width, height, depth, maxima);
	time(&end);
	std::cout << "Found " << maxima.size() << " nuclei." << std::endl;
	std::cout << "findMaxima took " << difftime(end, start) << " seconds." << std::endl;

	std::cout << "computing gradient field...";
	time(&start);
	gradientField3d(gaussian_result, width, height, depth, gradientField);
	std::cout << "done." << std::endl;
	time(&end);
	std::cout << "Gradient field took: " << difftime(end, start) << " seconds." << std::endl;

	std::cout << "segmenting blobs...";
	nuclei.resize(maxima.size());
	time(&start);
	concurrency::static_partitioner partitioner;

	concurrency::parallel_for(0, (int)maxima.size(), [&maxima, &nuclei, &width, &height, &depth, gradientField](int j) {
		Nucleus* blob = new Nucleus(1e7, 3e1);
		blob->id = j;
		blob->points.reserve(200000);
		blob->boundary.reserve(20000);
		blob->local_max = maxima.at(j);
		segment_blob(blob->points, blob->boundary, maxima.at(j), gradientField, width, height, depth);
		nuclei.at(j) = blob;
		}, partitioner);
	time(&end);

	std::cout << "done." << std::endl;
	std::cout << "blob segmentation took: " << difftime(end, start) << " seconds." << std::endl;
}


void findDots(uint16_t* stack,
	uint16_t* median,
	float* gaussian,
	std::vector<std::tuple<int, int, int>>& maxima,
	int threshold_lower,
	int threshold_higher,
	int sigmaxy,
	int sigmaz,
	int width, int height, int depth) {

	std::cout << "median...";
	medianFilter3x3(stack, width, height, depth, median);
	std::cout << "done." << std::endl;

	uint16_t* pointer = median;
	for (int z = 0; z < depth; z++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if ((*pointer) <= threshold_lower) {
					*pointer = 0;
				}
				else if ((*pointer) > threshold_higher) {
					*pointer = threshold_higher - threshold_lower;
				}
				else {
					*pointer = (*pointer) - threshold_lower;
				}
				pointer++;
			}
		}
	}
	std::cout << "computing gaussian...";
	gaussian_filter3D_parallel(median, width, height, depth, sigmaxy, sigmaz, gaussian);
	std::cout << "done." << std::endl;

	findMaxima(gaussian, width, height, depth, maxima);
	std::cout << "Found " << maxima.size() << " dots." << std::endl;
}


std::tuple<int, int> closestBlob(std::tuple<int, int, int> point, std::vector<Nucleus*> blobs, int celldiameter) {
	int n = 0;
	Blob* nuc = blobs.at(0);
	int min_dist = deltaXSq(point, nuc->points.at(0));
	for (int j = 1; j < nuc->points.size(); j++) {
		int dist = deltaXSq(point, nuc->points.at(j));
		if (dist < min_dist) {
			min_dist = dist;
		}
	}
	for (int k = 1; k < blobs.size(); k++) {
		nuc = blobs.at(k);
		// skip if we are more than 3 cell diameters away from cell center
		if (deltaXSq(point, nuc->local_max) > 9 * celldiameter * celldiameter) continue;
		for (int j = 0; j < nuc->points.size(); j++){
				int dist = deltaXSq(point, nuc->points.at(j));
				if (dist < min_dist) {
					min_dist = dist;
					n = k;
				}
			}
			if (min_dist == 0) break;
		}
		return std::make_tuple(min_dist, n);
}


cv::Mat closeUpCellDotMaxProjectImg(Nucleus* nuc, uint16_t* stack, uint16_t* median594, uint16_t* median647,
	float blue_pixel_lower,
	float blue_pixel_upper,
	float green_pixel_lower,
	float green_pixel_upper,
	float red_pixel_lower,
	float red_pixel_upper,
	float white_pixel_lower,
	float white_pixel_upper,
	int width, int height, int depth) {

	int nx, ny, nz, i;
	std::tie(nx, ny, nz) = nuc->local_max;
	i = nz;
	cv::Mat img(401, 401, CV_16UC3, cv::Scalar(0, 0, 0));

	// we are going to look at range from maxima-150 to maxima+150
	for (int rx = -200; rx <= 200; rx++) {
		if (nx + rx < 0 || nx + rx >= width) {
			continue;
		}
		for (int ry = -200; ry <= 200; ry++) {
			if (ny + ry < 0 || ny + ry >= height) {
				continue;
			}
			BGR& bgr = img.ptr<BGR>(ry + 200)[rx + 200];

			for (int zh = i - 50; zh < i + 50; zh++) {
				if (zh < 0) continue;
				if (zh >= depth) break;
				uint16_t green = median594[(nx + rx) + (ny + ry) * width + width * height * zh];
				if (bgr.green < 255 * 255 && bgr.green < (green - green_pixel_lower) * (255 * 255 / green_pixel_upper)) {
					if (green > green_pixel_upper) bgr.green = 255 * 255;
					else bgr.green = (green - green_pixel_lower) * (255 * 255 / green_pixel_upper);
				}
				uint16_t blue = stack[(nx + rx) + (ny + ry) * width + width * height * zh];

				if (bgr.blue < 255 * 255 && bgr.blue < (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper)) {
					if (blue > blue_pixel_upper) bgr.blue = 255 * 255;
					else bgr.blue = (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper);
				}

				uint16_t red = median647[(nx + rx) + (ny + ry) * width + width * height * zh];
				if (bgr.red < 255 * 255 && bgr.red < (red - red_pixel_lower) * (255 * 255 / red_pixel_upper)) {
					if (red > red_pixel_upper) bgr.red = 255 * 255;
					else bgr.red = (red - red_pixel_lower) * (255 * 255 / red_pixel_upper);
				}
			}
		}
	}

	// for all the dots, show their outline
	for (int j = 0; j < nuc->close_points594.size(); j++) {
		int dx, dy, dz;
		std::tie(dx, dy, dz) = nuc->close_points594.at(j);
		cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(0, 255 * 255, 0), 1, 8, 0);
		cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);

	}

	// for all the dots, show their outline
	for (int j = 0; j < nuc->close_points640.size(); j++) {
		int dx, dy, dz;
		std::tie(dx, dy, dz) = nuc->close_points640.at(j);
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
	return img;
}

void closeUpCellDotMaxProject(std::vector<Nucleus*> nuclei, uint16_t* stack, uint16_t* median594, uint16_t* median647,
	float blue_pixel_lower,
	float blue_pixel_upper,
	float green_pixel_lower,
	float green_pixel_upper,
	float red_pixel_lower,
	float red_pixel_upper,
	float white_pixel_lower,
	float white_pixel_upper,
	int width, int height, int depth) {
	int i, k = 0, key = 0;
	do {
		Nucleus* nuc = nuclei.at(k);
		cv::Mat img = closeUpCellDotMaxProjectImg(nuc, stack, median594, median647,
			blue_pixel_lower,
			blue_pixel_upper,
			green_pixel_lower,
			green_pixel_upper,
			red_pixel_lower,
			red_pixel_upper,
			white_pixel_lower,
			white_pixel_upper,
			width, height, depth);

		cv::imshow("Display window", img);

		key = cv::waitKey(0);
		if (key == 'a') {
			do {
				k = (k + 1) % nuclei.size();
			} while (!nuclei.at(k)->validSize());
		}
		else if (key == 'd') {
			do {
				k = (nuclei.size() + k - 1) % nuclei.size();
			} while (!nuclei.at(k)->validSize());
		}

	} while (key != 27);

}

void closeUpCellAndDots(std::vector<Nucleus*> nuclei, uint16_t* stack, uint16_t* median594, uint16_t* median647,
	float blue_pixel_lower,
	float blue_pixel_upper,
	float green_pixel_lower,
	float green_pixel_upper,
	float red_pixel_lower,
	float red_pixel_upper,
	float white_pixel_lower,
	float white_pixel_upper,
	int width, int height, int depth) {
	int k = 0, key = 0;
	int i = 0;
	do {
		int nx, ny, nz;
		Nucleus* nuc = nuclei.at(k);
		std::tie(nx, ny, nz) = nuc->local_max;
		cv::Mat img(401, 401, CV_16UC3, cv::Scalar(0, 0, 0));

		// we are going to look at range from maxima-150 to maxima+150
		for (int rx = -200; rx <= 200; rx++) {
			if (nx + rx < 0 || nx + rx >= width) {
				continue;
			}
			for (int ry = -200; ry <= 200; ry++) {
				if (ny + ry < 0 || ny + ry >= height) {
					continue;
				}
				BGR& bgr = img.ptr<BGR>(ry + 200)[rx + 200];

				uint16_t green = median594[(nx + rx) + (ny + ry) * width + width * height * i];
				if (green < green_pixel_lower) bgr.green = 0;
				else if (green > green_pixel_upper) bgr.green = 255 * 255;
				else bgr.green = (green - green_pixel_lower) * (255 * 255 / green_pixel_upper);

				uint16_t blue = stack[(nx + rx) + (ny + ry) * width + width * height * i];
				if (blue < blue_pixel_lower) bgr.blue = 0;
				else if (blue > blue_pixel_upper) bgr.blue = 255 * 255;
				else bgr.blue = (blue - blue_pixel_lower) * (255 * 255 / blue_pixel_upper);

				uint16_t red = median647[(nx + rx) + (ny + ry) * width + width * height * i];
				if (red < red_pixel_lower) bgr.red = 0;
				else if (red > red_pixel_upper) bgr.red = 255 * 255;
				else bgr.red = (red - red_pixel_lower) * (255 * 255 / red_pixel_upper);
			}
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

		// for all the dots, show their outline
		for (int j = 0; j < nuc->close_points594.size(); j++) {
			int dx, dy, dz;
			std::tie(dx, dy, dz) = nuc->close_points594.at(j);

			if (((dz - i) < 5 && i - dz < 5)) {
				// if we are close enough
				cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(0, 255 * 255, 0), 1, 8, 0);
				cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);

			}

		}

		// for all the dots, show their outline
		for (int j = 0; j < nuc->close_points640.size(); j++) {
			int dx, dy, dz;
			std::tie(dx, dy, dz) = nuc->close_points640.at(j);

			if (((dz - i) < 5 && i - dz < 5)) {
				//cv::rectangle(img, cv::Point(dx - 5, dy - 5), cv::Point(dx + 5, dy + 5), cv::Scalar(0, 0, 255), cv::FILLED, cv::LINE_8);
				cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(0, 0, 255 * 255), 1, 8, 0);
				cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(255 * 255, 255 * 255, 255 * 255), 1, 8, 0);
			}

		}
		cv::imshow("Display window", img);
		key = cv::waitKey(0);
		if (key == 's') {
			i = (i + 1) % depth;
		}
		else if (key == 'w') {
			i = (depth + i - 1) % depth;
		}
		else if (key == 'a') {
			do {
				k = (k + 1) % nuclei.size();
			} while (!nuclei.at(k)->validSize());
			i = std::get<2>(nuclei.at(k)->local_max);
		}
		else if (key == 'd') {
			do {
				k = (nuclei.size() + k - 1) % nuclei.size();
			} while (!nuclei.at(k)->validSize());
			i = std::get<2>(nuclei.at(k)->local_max);
		}
		else {
			//blink = !blink;
		}
	} while (key != 27);
}