#include <random>
#include <vector>
#include "pcg/pcg_basic.h"
#include <string>
#include "csv.h"
#include "mathutils.h"
#include "leastsquaresfit.h"

// TODO: temp
//#define DETERMINISTIC() false
#define DETERMINISTIC() true

// The size of the list of random numbers output
static const size_t c_numberCount = 10000000;

// TODO: uncomment the above
//static const size_t c_numberCount = 200;

// Bucket count of histogram that makes the PDF and CDF
static const size_t c_histogramBucketCount = 1024;

float PCGRandomFloat01(pcg32_random_t& rng)
{
	return ldexpf((float)pcg32_random_r(&rng), -32);
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
	}

	// TODO: polynomial fit the cdf and see how well if fits / does

	// Put the CDF into the CD Fcsv
	columnIndex = (int)CDFcsv.size();
	CDFcsv.resize(columnIndex + 4);
	CDFcsv[columnIndex].label = std::string(label) + " CDF";
	CDFcsv[columnIndex].values = CDF;

	// Fit the CDF with a polynomial
	LeastSquaresPolynomialFit<2> fit2;
	LeastSquaresPolynomialFit<3> fit3;
	LeastSquaresPolynomialFit<4> fit4;
	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float x = float(i) / float(CDF.size());
		float y = CDF[i];
		fit2.AddPoint(x, y);
		fit3.AddPoint(x, y);
		fit4.AddPoint(x, y);
	}

	fit2.CalculateCoefficients();
	fit3.CalculateCoefficients();
	fit4.CalculateCoefficients();

	CDFcsv[columnIndex + 1].label = std::string(label) + " Quadratic";
	CDFcsv[columnIndex + 2].label = std::string(label) + " Cubic";
	CDFcsv[columnIndex + 3].label = std::string(label) + " Quartic";

	CDFcsv[columnIndex + 1].values.resize(CDF.size());
	CDFcsv[columnIndex + 2].values.resize(CDF.size());
	CDFcsv[columnIndex + 3].values.resize(CDF.size());

	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float x = float(i) / float(CDF.size() - 1);
		CDFcsv[columnIndex + 1].values[i] = fit2.Evaluate(x);
		CDFcsv[columnIndex + 2].values[i] = fit3.Evaluate(x);
		CDFcsv[columnIndex + 3].values[i] = fit4.Evaluate(x);
	}

	// TODO: use the fits above to make the data uniform
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
