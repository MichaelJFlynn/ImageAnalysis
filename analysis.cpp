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
// the main hpp file
#include "ImageAnalysis.hpp"



int main()
{
	// files
	std::array<char*, 4> construct_dirs = {
		//"MouseBrain2021-10-21\\0x\\0x_slice1_",
		"MouseBrain2021-10-21\\1x\\1x_",
		//"MouseBrain2021-10-21\\10x\\10x_",
		//"MouseBrain2021-10-21\\endogenous_0x\\endogenous_0x_"
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
		int pos = 8;
		
		char dapi_filename[500];
		sprintf(dapi_filename, "%s#%d_New.tif", dir, pos);

		FILE* dapi_file;
		if((dapi_file = fopen(dapi_filename, "r")) && pos < 10) {

		//while ((dapi_file = fopen(dapi_filename, "r")) && pos < 10) {
			fclose(dapi_file);
			std::cout << "Working on " << dapi_filename << std::endl;

			// fill out filenames
			char file_488[500];
			char file_594[500];
			char file_647[500];
			sprintf(file_488, "%s#%d_New_1.tif", dir, pos);
			sprintf(file_594, "%s#%d_New_2.tif", dir, pos);
			sprintf(file_647, "%s#%d_New_4.tif", dir, pos);


			int threshold_405_lower = 200;
			int threshold_405_higher = 390;
			//std::cout << "please observe this dog..." << std::endl;
			std::cout << "Loading DAPI...";
			//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
			uint16_t* stack = new uint16_t[2048 * 2048 * 201];
			loadTiff(stack, dapi_filename, 2048, 2048, 201);
			std::cout << " done." << std::endl;
			time_t start, end;

			uint16_t* filtered = new uint16_t[2048 * 2048 * 201];
			std::cout << "Computing median filter...";
			time(&start);
			medianFilter3x3(stack, 2048, 2048, 201, filtered);
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
			gaussian_filter3D_parallel(filtered, 2048, 2048, 201, sigma, 10, gaussian_result);
			std::cout << "done." << std::endl;

			time(&end);
			std::cout << "gaussian filter took " << difftime(end, start) << " seconds." << std::endl;

			std::vector<std::tuple<int, int, int>> maxima;
			time(&start);
			findMaxima(gaussian_result, 2048, 2048, 201, maxima);
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
			gradientField3d(gaussian_result, 2048, 2048, 201, gradientField);
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
				Nucleus* blob = new Nucleus(1e6, 1e5);
				blob->id = j;
				blob->points.reserve(200000);
				blob->boundary.reserve(20000);
				blob->local_max = maxima.at(j);
				segment_blob(blob->points, blob->boundary, maxima.at(j), gradientField, 2048, 2048, 201);
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
			loadTiff(stack488, file_488, 2048, 2048, 201);
			std::cout << "done." << std::endl;

			std::cout << "median 488...";
			uint16_t* median488 = new uint16_t[2048 * 2048 * 201];
			medianFilter3x3(stack488, 2048, 2048, 201, median488);
			std::cout << "done." << std::endl;
			delete[] stack488;

			
			pointer = median488;
			for (int z = 0; z < 201; z++) {
				for (int y = 0; y < 2048; y++) {
					for (int x = 0; x < 2048; x++) {
						if ((*pointer) <= 130) {
							*pointer = 0;
						}
						else {
							*pointer = (*pointer) - 130; // this should be 160 to accurately determine GFP signal over noise in nucleii.
						}
						pointer++;
					}
				}
			} // zero it out

			// compute gaussian
			std::cout << "computing 488 gaussian...";
			float* gaussian_488 = new float[2048 * 2048 * 201];
			gaussian_filter3D_parallel(median488, 2048, 2048, 201, 2, 1, gaussian_488);
			std::cout << "done." << std::endl;

			// find maxima
			std::vector<std::tuple<int, int, int>> maxima488;
			findMaxima(gaussian_488, 2048, 2048, 201, maxima488);
			std::cout << "Found " << maxima488.size() << " dots." << std::endl;



			// get gradient field
			// not doing this for now, not really interested in segmenting by size


			delete[] gaussian_488;
			// segment dots

			//std::cout << "median filter 488...";
			//uint16_t* median488 = new uint16_t[2048 * 2048 * 201];
			//medianFilter3x3(stack488, median488);
			//std::cout << "done." << std::endl;


			//std::cout << "please observe this dog..." << std::endl;
			std::cout << "Loading 594...";
			//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
			uint16_t* stack594 = new uint16_t[2048 * 2048 * 201];
			loadTiff(stack594, file_594, 2048, 2048, 201);
			std::cout << " done." << std::endl;

			std::cout << "median 594...";
			uint16_t* median594 = new uint16_t[2048 * 2048 * 201];
			medianFilter3x3(stack594, 2048, 2048, 201, median594);
			std::cout << "done." << std::endl;
			delete[] stack594;

			pointer = median594;
			for (int z = 0; z < 201; z++) {
				for (int y = 0; y < 2048; y++) {
					for (int x = 0; x < 2048; x++) {
						if ((*pointer) <= 120) {
							*pointer = 0;
						}
						else {
							*pointer = (*pointer) - 120;
						}
						pointer++;
					}
				}
			}
			std::cout << "computing 594 gaussian...";
			float* gaussian_594 = new float[2048 * 2048 * 201];
			gaussian_filter3D_parallel(median594, 2048, 2048, 201, 2, 1, gaussian_594);
			std::cout << "done." << std::endl;

			std::vector<std::tuple<int, int, int>> maxima594;
			findMaxima(gaussian_594, 2048, 2048, 201, maxima594);
			std::cout << "Found " << maxima594.size() << " dots." << std::endl;

			std::cout << "computing 594 gradient field...";
			time(&start);
			gradientField3d(gaussian_594, 2048, 2048, 201, gradientField);
			std::cout << "done." << std::endl;
			time(&end);
			std::cout << "Gradient field took: " << difftime(end, start) << " seconds." << std::endl;
			delete[] gaussian_594;

			std::cout << "segmenting 594 dots...";
			std::vector<Dot*> dots594;
			dots594.resize(maxima594.size());
			concurrency::parallel_for(0, (int)maxima594.size(), [&maxima594, &gradientField, &dots594](int j) {
				Dot* blob = new Dot(5000, 1);
				blob->id = j;
				blob->points.reserve(2000);
				blob->boundary.reserve(200);
				blob->local_max = maxima594.at(j);
				segment_blob(blob->points, blob->boundary, maxima594.at(j), gradientField, 2048, 2048, 201);
				dots594.at(j) = blob;
				}, partitioner);
			std::cout << "done." << std::endl;
			std::cout << "Segmented " << dots594.size() << " dots." << std::endl;

			//std::cout << "please observe this dog..." << std::endl;
			std::cout << "Loading 640...";
			//std::vector<cv::Mat*>* mats = load_tiff("endogenous_0x_#1_New.tif");
			uint16_t* stack647 = new uint16_t[2048 * 2048 * 201];
			loadTiff(stack647, file_647, 2048, 2048, 201);
			std::cout << " done." << std::endl;

			std::cout << "median 640...";
			uint16_t* median640 = new uint16_t[2048 * 2048 * 201];
			medianFilter3x3(stack647, 2048, 2048, 201, median640);
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
			gaussian_filter3D_parallel(median640, 2048, 2048, 201, 2, 1, gaussian_647);
			std::cout << "done." << std::endl;

			std::vector<std::tuple<int, int, int>> maxima647;
			findMaxima(gaussian_647, 2048, 2048, 201, maxima647);
			std::cout << "Found " << maxima647.size() << " dots." << std::endl;

			std::cout << "computing 647 gradient field...";
			time(&start);
			gradientField3d(gaussian_647, 2048, 2048, 201, gradientField);
			std::cout << "done." << std::endl;
			time(&end);
			std::cout << "Gradient field took: " << difftime(end, start) << " seconds." << std::endl;
			delete[] gaussian_647;

			std::cout << "segmenting 640 dots...";
			std::vector<Dot*> dots647;
			dots647.resize(maxima647.size());
			concurrency::parallel_for(0, (int)maxima647.size(), [&maxima647, &gradientField, &dots647](int j) {
				Dot* blob = new Dot(5000, 1);
				blob->id = j;
				blob->points.reserve(2000);
				blob->boundary.reserve(200);
				blob->local_max = maxima647.at(j);
				segment_blob(blob->points, blob->boundary, maxima647.at(j), gradientField, 2048, 2048, 201);
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
					sumgreen += median488[xg + 2048 * yg + 2048 * 2048 * zg];
				}
				fprintf(nucleii_csv, "%d,%d,%d,%d,%d,%s,%d,%s,%d\n", nuc->id, max_x, max_y, max_z, nuc->points.size(), nuc->validSize() ? "TRUE" : "FALSE", sumgreen, dir, pos);
			}
			fclose(nucleii_csv);
			
			// make the dots csv
			// obviously we need x,y,z,channel, 
			// also closest nuclei id, and x,y,z
			int id = 0;
			dots_csv = fopen("dots.csv", "a");
			for (int i = 0; i < maxima488.size(); i++) {
				int x, y, z;
				std::tie(x, y, z) = maxima488.at(i);

				fprintf(dots_csv, "%d,%d,%d,%d,,,,,488,%s,%d\n", id, x, y, z, dir, pos);
				id++;
			}


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
					fprintf(dots_csv, "%d,%d,%d,%d,%d,%d,%d,%d,594,%s,%d\n", id, x, y, z, cn_id, cn_x, cn_y, cn_z, dir, pos);
					id++;
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
					fprintf(dots_csv, "%d,%d,%d,%d,%d,%d,%d,%d,640,%s,%d\n", id, x, y, z, cn_id, cn_x, cn_y, cn_z, dir, pos);
					id++;
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
			delete[] median488;
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