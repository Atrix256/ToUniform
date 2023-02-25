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
	std::vector<float> values;
};
typedef std::vector<Column> CSV;

float PCGRandomFloat01(pcg32_random_t& rng)
{
	return ldexpf((float)pcg32_random_r(&rng), -32);
}

std::vector<float> Convolve(const std::vector<float>& A, const std::vector<float>& B)
{
	const int sizeA = int(A.size());
	const int sizeB = int(B.size());
	const int sizeOut = sizeA + sizeB - 1;

	std::vector<float> out(sizeOut, 0.0f);

	for (int outIndex = 0; outIndex < sizeOut; ++outIndex)
	{
		int indexA = outIndex - sizeB + 1;
		int indexB = sizeB - 1;

		if (indexA < 0)
		{
			indexB += indexA;
			indexA = 0;
		}

		while (indexA < sizeA && indexB >= 0)
		{
			out[outIndex] += A[indexA] * B[indexB];
			indexA++;
			indexB--;
		}
	}

	return out;
}

void KernelTest(const char* label, pcg32_random_t& rng, CSV& csv, const std::vector<float>& kernel)
{
	printf("%s\n", label);

	// reserve space in the CSV for this data
	int columnIndex = (int)csv.size();
	csv.resize(columnIndex + 2);
	csv[columnIndex].label = label;
	csv[columnIndex + 1].label = std::string(label) + "_ToUniform";

	// make white noise
	printf("  Generating\n");
	std::vector<float> whiteNoise(c_numberCount);
	for (float& f : whiteNoise)
		f = PCGRandomFloat01(rng);

	// filter the white noise, truncate it, and put it into the csv
	printf("  Filterting\n");
	csv[columnIndex].values = Convolve(whiteNoise, kernel);
	csv[columnIndex].values.resize(c_numberCount);

	// TODO: make PDF and CDF and then put the filtered noise through that.
	csv[columnIndex + 1].values = csv[columnIndex].values;
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

	// make noise and add to csv
	CSV csv;
	KernelTest("Box2RedNoise", rng, csv, { 1.0f, 1.0f });
	KernelTest("Box2BlueNoise", rng, csv, { 1.0f, -1.0f });

	// TODO: more types of noise, and gaussian filtered.

	// write to a csv so we can make histograms
	printf("\nWriting CSV...\n");
	FILE* file = nullptr;
	fopen_s(&file, "out.csv", "wb");

	// write header
	bool first = true;
	for (const Column& column : csv)
	{
		fprintf(file, "%s\"%s\"", first ? "" : ",", column.label.c_str());
		first = false;
	}
	fprintf(file, "\n");

	// write data
	for (size_t i = 0; i < c_numberCount; ++i)
	{
		bool first = true;
		for (const Column& column : csv)
		{
			fprintf(file, "%s\"%f\"", first ? "" : ",", column.values[i]);
			first = false;
		}
		fprintf(file, "\n");
	}

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
- since it's now a noise type and then ToUniform, we should make the output have 2 columns, and N rows where N is the number of noise types we have.
 - could also do 4, 6 or 8 columns instead to make the image more square

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
