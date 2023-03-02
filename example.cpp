#include "ImageAnalysis.hpp"


int main()
{
	const char* oligoDT_filename = "example_405.tif";
	std::cout << "Checking if " << oligoDT_filename << " exists." << std::endl;
	int width, height, depth;
	getTiffDimensions(oligoDT_filename, width, height, depth);
	printf("Image width: %d, height: %d, depth: %d\n", width, height, depth);

	// load 
	std::cout << "Loading 405...";
	uint16_t* stack = new uint16_t[width * height * depth];
	loadTiff(stack, oligoDT_filename,width,height, depth);
	std::cout << " done." << std::endl;

	//scanInt(stack, width, height, depth);

	uint16_t* filtered = new uint16_t[width * height * depth];
	float* gaussian_result = new float[width * height * depth];
	float3* gradientField = new float3[width * height * depth];
	std::vector<Nucleus*> nuclei;
	std::vector<std::tuple<int, int, int>> maxima;
	int threshold_405_higher = 150;
	int threshold_405_lower = 110;
	int sigmaxy = 20;
	int sigmaz = 10;
	segmentNuclei(stack,
		filtered,
		gaussian_result,
		gradientField,
		nuclei,
		maxima,
		threshold_405_lower,
		threshold_405_higher,
		sigmaxy,
		sigmaz,
		width, height, depth);

	int num_good = 0;
	for (int i = 0; i < nuclei.size(); i++) {
		if (nuclei.at(i)->validSize()) {
			num_good++;
		}
	}
	std::cout << "Found " << num_good << " out of " << nuclei.size() << " maxima to be good..." << std::endl;

	// find 594 dots
	uint16_t* stack594 = new uint16_t[width * height * depth];
	uint16_t* median594 = new uint16_t[width * height * depth];
	float* gaussian594 = new float[width * height * depth];
	std::vector<std::tuple<int, int, int>> maxima594;
	int threshold_594_lower = 200;
	int threshold_594_higher = 1e6;
	int sigmaxy_594 = 2;
	int sigmaz_594 = 1;

	std::cout << "Loading 594...";
	char* file_594 = "example_594.tif";
	loadTiff(stack594, file_594, width, height, depth);
	std::cout << " done." << std::endl;

	std::cout << "Finding 594 dots..." << std::endl;

	findDots(stack594,
		median594,
		gaussian594,
		maxima594,
		threshold_594_lower,
		threshold_594_higher,
		sigmaxy_594,
		sigmaz_594,
		width, height, depth);

	// 640 section
	uint16_t* stack647 = new uint16_t[width * height * depth];
	uint16_t* median647 = new uint16_t[width * height * depth];
	float* gaussian647 = new float[width * height * depth];
	std::vector<std::tuple<int, int, int>> maxima647;
	int threshold_647_lower = 200;
	int threshold_647_higher = 1e6;
	int sigmaxy_647 = 2;
	int sigmaz_647 = 1;

	std::cout << "Loading 640...";
	char* file_647 = "example_640.tif";
	loadTiff(stack647, file_647, width, height, depth);
	std::cout << " done." << std::endl;


	std::cout << "Finding 640 dots..." << std::endl;

	findDots(stack647,
		median647,
		gaussian647,
		maxima647,
		threshold_647_lower,
		threshold_647_higher,
		sigmaxy_647,
		sigmaz_647,
		width, height, depth);

	// Match dots to nuclei
	std::cout << "Matching dots to cells...";
	int threshold_distSq = 10 * 10;
	int cell_diameter = 100;
	// trying to match up dots with nuclei.
	for (int i = 0; i < maxima647.size(); i++) {
		int min_dist, n;
		std::tie(min_dist, n) = closestBlob(maxima647.at(i), nuclei, cell_diameter);
		if (min_dist < threshold_distSq) {
			nuclei.at(n)->close_points640.push_back(maxima647.at(i));
		}
	}

	for (int i = 0; i < maxima594.size(); i++) {
		int min_dist, n;
		std::tie(min_dist, n) = closestBlob(maxima594.at(i), nuclei, cell_diameter);
		if (min_dist < threshold_distSq) {
			nuclei.at(n)->close_points594.push_back(maxima594.at(i));
		}
	}
	std::cout << "done." << std::endl;

	// We can display a bunch of things now 
	const float blue_pixel_lower = 110;
	const float blue_pixel_upper = 200;
	const float green_pixel_lower = 0;
	const float green_pixel_upper = 100;
	const float red_pixel_lower = 0;
	const float red_pixel_upper = 100;
	const float white_pixel_lower = 0;
	const float white_pixel_upper = 1;

	scanInt(stack, width, height, depth);
	scanInt(filtered, width, height, depth);
	

	displayBlobsInt(stack, width, height, depth, nuclei);

	closeUpCellAndDots(nuclei, stack, median594, median647,
		blue_pixel_lower,
		blue_pixel_upper,
		green_pixel_lower,
		green_pixel_upper,
		red_pixel_lower,
		red_pixel_upper,
		white_pixel_lower,
		white_pixel_upper,
		width, height, depth);

	closeUpCellDotMaxProject(nuclei, stack, median594, median647,
		blue_pixel_lower,
		blue_pixel_upper,
		green_pixel_lower,
		green_pixel_upper,
		red_pixel_lower,
		red_pixel_upper,
		white_pixel_lower,
		white_pixel_upper,
		width, height, depth);

	std::cout << "Loading 488...";
	char* file_488 = "example_488.tif";
	uint16_t* stack488 = new uint16_t[width * height * depth];
	loadTiff(stack488, file_488, width, height, depth);
	std::cout << " done." << std::endl;

	std::cout << "Loading 561...";
	char* file_561 = "example_561.tif";
	uint16_t* stack561 = new uint16_t[width * height * depth];
	loadTiff(stack561, file_561, width, height, depth);
	std::cout << " done." << std::endl;

	std::cout << "Writing data to csv files...";

	// write dot data to file
	FILE* dots_csv = fopen("dots_example.csv", "w");
	FILE* nuclei_csv = fopen("nuclei_example.csv", "w");
	fprintf(dots_csv, "x,y,z,cn_id,cn_x,cn_y,cn_z,channel\n");
	fprintf(nuclei_csv, "cn_id,cn_x,cn_y,cn_z,size,validSize,sum488,sum561\n");
	for (int i = 0; i < nuclei.size(); i++) {
		int cn_x, cn_y, cn_z;
		Nucleus nuc = *nuclei.at(i);
		std::tie(cn_x, cn_y, cn_z) = nuc.local_max;
		int sum488 = 0;
		int sum561 = 0;
		for (int j = 0; j < nuc.points.size(); j++) {
			int x, y, z;
			std::tie(x, y, z) = nuc.points.at(j);
			sum488 += stack488[x + width * y + width * height * z];
			sum561 += stack561[x + width * y + width * height * z];
		}

		fprintf(nuclei_csv, "%d,%d,%d,%d,%d,%d,%d,%d\n", i, cn_x, cn_y, cn_z, nuc.points.size(), nuc.validSize(), sum488, sum561);


		for (int j = 0; j < nuc.close_points594.size(); j++) {
			int x, y, z;
			std::tie(x, y, z) = nuc.close_points594.at(j);
			fprintf(dots_csv, "%d,%d,%d,%d,%d,%d,%d,%d\n", x, y, z, i, cn_x, cn_y, cn_z, 594);

		}
		for (int j = 0; j < nuc.close_points640.size(); j++) {
			int x, y, z;
			std::tie(x, y, z) = nuc.close_points640.at(j);
			fprintf(dots_csv, "%d,%d,%d,%d,%d,%d,%d,%d\n", x, y, z, i, cn_x, cn_y, cn_z, 640);
		}

	}
	fclose(dots_csv);
	fclose(nuclei_csv);
	std::cout << "Done." << std::endl;



	// if in a loop, make sure to delete allocated memory (everying line with 'new')
	for (int i = 0; i < nuclei.size(); i++) {
		delete nuclei.at(i);
	}
	delete[] stack;
	delete[] filtered;
	delete[] gaussian_result;
	delete[] gradientField;
	delete[] stack594;
	delete[] median594;
	delete[] gaussian594;
	delete[] stack647;
	delete[] median647;
	delete[] gaussian647;
	delete[] stack488;
	delete[] stack561;
}