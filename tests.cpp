#include "ImageAnalysis.hpp"

typedef std::tuple<int, int, int> coord3d;

class Test {
public: 
	std::string title;
	bool failed;
	Test(std::string t) : title(t), failed(false) {
	}
};

int main() {
	std::cout << "Testing Image Analysis..." << std::endl;
	std::vector<Test> tests;

	// test gaussian filter
	Test gaussianTest("Testing if gaussian filter produces accurate results on delta functions.");
	std::cout << gaussianTest.title << "..";
	uint16_t* test_stack = new uint16_t[100 * 100 * 100];
	std::fill(test_stack, test_stack + 100*100*100, 0);
	coord3d delta_coordinates[4] = {
		std::make_tuple(50, 50, 50),
		std::make_tuple(40, 50, 50),
		std::make_tuple(50, 40, 50),
		std::make_tuple(50, 50, 40),

	};
	for(int i = 0; i < 4; i++) {
		int x, y, z;
		std::tie(x, y, z) = delta_coordinates[i];
		test_stack[x + 100 * y + 100*100*z] = 1;

	}
	float* result = new float[100 * 100 * 100];
	gaussian_filter3D_parallel(test_stack, 100, 100, 100, 10, 10, result);
	float norm = 1.0 / sqrt(2 * M_PI * 10 * 10);

	bool success = true;
	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			for (int k = 0; k < 100; k++) { 
				float expected_value = 0;
				for (int l = 0; l < 4; l++) {
					int x, y, z;
					std::tie(x, y, z) = delta_coordinates[l];
					float dx = x - i;
					float dy = y - j;
					float dz = z - k;
					expected_value += exp(-(dx * dx + dy * dy + dz * dz) / (2 * 10 * 10)) * norm * norm * norm;
				}
				
				if (abs(result[i + 100 * j + 100 * 100 * k] - expected_value) > norm*norm*norm* 1e-6) {
					gaussianTest.failed = true;
				}
			}
		}
	}
	//std::cout << norm * norm * norm << " " << result[50 + 100 * 50 + 100 * 100 * 50];
	if (gaussianTest.failed) {
		std::cout << "FAILED." << std::endl;
	}
	else {
		std::cout << "passed." << std::endl;

	}
	tests.push_back(gaussianTest);




	// test resizing.


	std::cout << "Tests done." << std::endl;
	return 0;
}