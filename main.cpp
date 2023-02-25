#include <random>
#include <vector>
#include "pcg/pcg_basic.h"
#include <string>

#define DETERMINISTIC() false

// TODO: back to this
static const size_t c_numberCount = 10000000;
//static const size_t c_numberCount = 200;

struct Column
{
	std::string label;
	std::vector<float> value;
};
typedef std::vector<Column> CSV;

float PCGRandomFloat01(pcg32_random_t& rng)
{
	return ldexpf((float)pcg32_random_r(&rng), -32);
}

template <typename LAMBDA>
void KernelTest(const char* label, CSV& csv, const LAMBDA& lambda)
{
	int columnIndex = (int)csv.size();
	csv.resize(columnIndex + 2);
	csv[columnIndex].label = label;
	csv[columnIndex + 1].label = std::string(label) + "ToUniform";
}

int main(int argc, char** argv)
{
	pcg32_random_t rng;
#if !DETERMINISTIC()
	std::random_device rd;
	pcg32_srandom_r(&rng, rd(), 0);
#else
	pcg32_srandom_r(&rng, 0xa000b800, 0);
#endif

	// make triangular red noise through addition (low pass filter)
	printf("order2RedNoise\n");
	std::vector<float> order2RedNoise(c_numberCount);
	{
		float a[2] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 2] = PCGRandomFloat01(rng);
			order2RedNoise[i] = a[0] + a[1];
			order2RedNoise[i] = order2RedNoise[i] / 2.0f;
		}
	}

	// make triangular blue noise through subtraction (high pass filter)
	printf("order2BlueNoise\n");
	std::vector<float> order2BlueNoise(c_numberCount);
	{
		float a[2] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 2] = PCGRandomFloat01(rng);
			order2BlueNoise[i] = a[i % 2] - a[(i + 1) % 2];
			order2BlueNoise[i] = (1.0f + order2BlueNoise[i]) / 2.0f;
		}
	}

	// order 3 blue noise through subtraction (high pass filter)
	printf("order3BlueNoise\n");
	std::vector<float> order3BlueNoise(c_numberCount);
	{
		float a[3] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 3] = PCGRandomFloat01(rng);

			bool plus = true;
			order3BlueNoise[i] = 0.0f;
			for (int j = 0; j < 3; ++j)
			{
				order3BlueNoise[i] += a[(i + j) % 3] * (plus ? 1.0f : -1.0f);
				plus = !plus;
			}

			order3BlueNoise[i] = (order3BlueNoise[i] + 1.0f) / 3.0f;
		}
	}

	// order 4 blue noise through subtraction (high pass filter)
	printf("order4BlueNoise\n");
	std::vector<float> order4BlueNoise(c_numberCount);
	{
		float a[4] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 4] = PCGRandomFloat01(rng);

			bool plus = true;
			order4BlueNoise[i] = 0.0f;
			for (int j = 0; j < 4; ++j)
			{
				order4BlueNoise[i] += a[(i + j) % 4] * (plus ? 1.0f : -1.0f);
				plus = !plus;
			}

			order4BlueNoise[i] = (order4BlueNoise[i] + 2.0f) / 4.0f;
		}
	}

	// order 5 blue noise through subtraction (high pass filter)
	printf("order5BlueNoise\n");
	std::vector<float> order5BlueNoise(c_numberCount);
	{
		float a[5] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 5] = PCGRandomFloat01(rng);

			bool plus = true;
			order5BlueNoise[i] = 0.0f;
			for (int j = 0; j < 5; ++j)
			{
				order5BlueNoise[i] += a[(i + j) % 5] * (plus ? 1.0f : -1.0f);
				plus = !plus;
			}

			order5BlueNoise[i] = (order5BlueNoise[i] + 2.0f) / 5.0f;
		}
	}

	// order 3 red noise through addition (low pass filter)
	printf("order3RedNoise\n");
	std::vector<float> order3RedNoise(c_numberCount);
	{
		float a[3] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 3] = PCGRandomFloat01(rng);
			order3RedNoise[i] = a[0] + a[1] + a[2];
			order3RedNoise[i] = order3RedNoise[i] / 3.0f;
		}
	}

	// order 4 red noise through addition (low pass filter)
	printf("order4RedNoise\n");
	std::vector<float> order4RedNoise(c_numberCount);
	{
		float a[4] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 4] = PCGRandomFloat01(rng);
			order4RedNoise[i] = a[0] + a[1] + a[2] + a[3];
			order4RedNoise[i] = order4RedNoise[i] / 4.0f;
		}
	}

	// order 5 red noise through addition (low pass filter)
	printf("order5RedNoise\n");
	std::vector<float> order5RedNoise(c_numberCount);
	{
		float a[5] = { PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng), PCGRandomFloat01(rng) };
		for (size_t i = 0; i < c_numberCount; ++i)
		{
			a[i % 5] = PCGRandomFloat01(rng);
			order5RedNoise[i] = a[0] + a[1] + a[2] + a[3] + a[4];
			order5RedNoise[i] = order5RedNoise[i] / 5.0f;
		}
	}

	// write to a csv so we can make histograms
	printf("\nWriting CSV...\n");
	FILE* file = nullptr;
	fopen_s(&file, "out.csv", "wb");
	fprintf(file,
		"\"order2RedNoise\",\"Order3RedNoise\",\"Order4RedNoise\",\"Order5RedNoise\","
		"\"order2BlueNoise\",\"Order3BlueNoise\",\"Order4BlueNoise\",\"Order5BlueNoise\""
		"\n"
	);
	for (size_t i = 0; i < c_numberCount; ++i)
		fprintf(file,
			"\"%f\",\"%f\",\"%f\",\"%f\","
			"\"%f\",\"%f\",\"%f\",\"%f\""
			"\n",
			order2RedNoise[i], order3RedNoise[i], order4RedNoise[i], order5RedNoise[i],
			order2BlueNoise[i], order3BlueNoise[i], order4BlueNoise[i], order5BlueNoise[i]
	);
	fclose(file);

	return 0;
}

/*
TODO:

- make a thing where you give a kernel and it samples to make PDF, then integrates to get CDF, and uses that to make uniform again.
 - CDF be stored as a LUT and also as least squared polynomial?
 - show resulting histogram and DFT
- better DFTs by averaging a bunch of smaller sections.
- remove the code in this file you don't need anymore

LANDFILL:

- irwin hall distribution CDF to convert back to uniform
 - LUT for irwin hall CDF, or least squares polynomial fit?
 - or turn the special cases of PDF piecewise polynomial into CDF
 - irwin hall cdf doesn't help weighted samples.

? what is irwin hall for weighted values? or aka uniform that is within different ranges
 * two make a trapezoid.
 * more make something more generalized.
 * maybe make a PDF through sampling, then integrate


*/
