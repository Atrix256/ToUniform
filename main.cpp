#include <random>
#include <vector>
#include "pcg/pcg_basic.h"
#include <string>

#define DETERMINISTIC() false

static const size_t c_numberCount = 10000000;

// Bucket count of histogram that makes the PDF and CDF
static const size_t c_histogramBucketCount = 1024;

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

float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
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

void KernelTest(const char* label, pcg32_random_t& rng, CSV& csv, CSV& CDFcsv, const std::vector<float>& kernel)
{
	printf("%s\n", label);

	// reserve space in the CSV for this data
	int columnIndex = (int)csv.size();
	csv.resize(columnIndex + 2);
	csv[columnIndex].label = label;
	csv[columnIndex + 1].label = std::string(label) + "_ToUniform";

	// make white noise
	std::vector<float> whiteNoise(c_numberCount);
	for (float& f : whiteNoise)
		f = PCGRandomFloat01(rng);

	// filter the white noise, truncate it, normalize it to [0,1] and put it into the csv
	csv[columnIndex].values = Convolve(whiteNoise, kernel);
	csv[columnIndex].values.resize(c_numberCount);
	float themin = csv[columnIndex].values[0];
	float themax = csv[columnIndex].values[0];
	for (float f : csv[columnIndex].values)
	{
		themin = std::min(themin, f);
		themax = std::max(themax, f);
	}
	for (float& f : csv[columnIndex].values)
		f = (f - themin) / (themax - themin);

	// Make a histogram
	std::vector<int> counts(c_histogramBucketCount, 0);
	for (float f : csv[columnIndex].values)
	{
		int bucket = std::min(int(f * float(c_histogramBucketCount)), (int)c_histogramBucketCount - 1);
		counts[bucket]++;
	}

	// Make a PDF and a CDF
	std::vector<float> PDF(c_histogramBucketCount);
	std::vector<float> CDF(c_histogramBucketCount);
	float lastCDFValue = 0.0f;
	for (int index = 0; index < c_histogramBucketCount; ++index)
	{
		PDF[index] = float(counts[index]) / float(c_numberCount);
		CDF[index] = lastCDFValue + PDF[index];
		lastCDFValue = CDF[index];
	}
	CDF.insert(CDF.begin(), 0.0f);
	CDF.push_back(1.0f);

	// Put the values through the CDF (inverted, inverted CDF) to make them be a uniform distribution
	csv[columnIndex + 1].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[columnIndex].values[index];
		float xindexf = std::min(x * float(c_histogramBucketCount), (float)(c_histogramBucketCount - 1));
		int xindex1 = int(xindexf);
		int xindex2 = std::min(xindex1 + 1, (int)c_histogramBucketCount - 1);
		float xindexfract = xindexf - std::floor(xindexf);

		float y1 = CDF[xindex1];
		float y2 = CDF[xindex2];

		csv[columnIndex + 1].values[index] = Lerp(y1, y2, xindexfract);

		// TODO: with low bucket counts, the zero bucket is low, and the last bucket is high. may be off by 0.5 bucket issue.
	}

	// TODO: after using fraction to lerp between buckets, may be able to drop bucket count.
	// TODO: first CDF value isn't 0. last is 1 because we force it to be. is this is a problem?
	// TODO: polynomial fit the cdf and see how well if fits / does

	// Put the CDF into the CD Fcsv
	columnIndex = (int)CDFcsv.size();
	CDFcsv.resize(columnIndex + 1);
	CDFcsv[columnIndex].label = std::string(label) + " CDF";
	CDFcsv[columnIndex].values = CDF;
}

void WriteCSV(const CSV& csv, const char* fileName)
{
	FILE* file = nullptr;
	fopen_s(&file, fileName, "wb");

	if (csv.size() > 0)
	{
		// write header
		size_t maxColSize = 0;
		bool first = true;
		for (const Column& column : csv)
		{
			fprintf(file, "%s\"%s\"", first ? "" : ",", column.label.c_str());
			maxColSize = std::max(maxColSize, column.values.size());
			first = false;
		}
		fprintf(file, "\n");

		// write data
		for (size_t i = 0; i < maxColSize; ++i)
		{
			bool first = true;
			for (const Column& column : csv)
			{
				if (column.values.size() > i)
					fprintf(file, "%s\"%f\"", first ? "" : ",", column.values[i]);
				else
					fprintf(file, "%s\"\"", first ? "" : ",");
				first = false;
			}
			fprintf(file, "\n");
		}
	}

	fclose(file);
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
	CSV csv, CDFcsv;
	KernelTest("Box2RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f });
	KernelTest("Box3RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f });
	KernelTest("Box4RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f, 1.0f });
	KernelTest("Box5RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f });
	KernelTest("Box2BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f });
	KernelTest("Box3BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f, 1.0f });
	KernelTest("Box4BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f, 1.0f, -1.0f });
	KernelTest("Box5BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f, 1.0f, -1.0f, 1.0f });

	// TODO: more types of noise, and gaussian filtered.

	printf("\nWriting CSVs...\n");
	WriteCSV(csv, "out.csv");
	WriteCSV(CDFcsv, "cdf.csv");

	return 0;
}

/*
TODO:

- better DFTs by averaging a bunch of smaller sections.

! may want to hold off on this post til you put out the paper?

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
