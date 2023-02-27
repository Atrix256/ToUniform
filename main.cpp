#include <random>
#include <vector>
#include "pcg/pcg_basic.h"
#include <string>
#include "csv.h"
#include "mathutils.h"
#include "leastsquaresfit.h"

#define DETERMINISTIC() false

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
	int csvcolumnIndex = (int)csv.size();
	csv.resize(csvcolumnIndex + 3);
	csv[csvcolumnIndex].label = label;
	csv[csvcolumnIndex + 1].label = std::string(label) + "_ToUniform";

	// make white noise
	std::vector<float> whiteNoise(c_numberCount);
	for (float& f : whiteNoise)
		f = PCGRandomFloat01(rng);

	// filter the white noise, truncate it, normalize it to [0,1] and put it into the csv
	csv[csvcolumnIndex].values = Convolve(whiteNoise, kernel);
	csv[csvcolumnIndex].values.resize(c_numberCount);
	float themin = csv[csvcolumnIndex].values[0];
	float themax = csv[csvcolumnIndex].values[0];
	for (float f : csv[csvcolumnIndex].values)
	{
		themin = std::min(themin, f);
		themax = std::max(themax, f);
	}
	for (float& f : csv[csvcolumnIndex].values)
		f = (f - themin) / (themax - themin);

	// Make a histogram
	std::vector<int> counts(c_histogramBucketCount, 0);
	for (float f : csv[csvcolumnIndex].values)
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
	csv[csvcolumnIndex + 1].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[csvcolumnIndex].values[index];
		float xindexf = std::min(x * float(c_histogramBucketCount), (float)(c_histogramBucketCount - 1));
		int xindex1 = int(xindexf);
		int xindex2 = std::min(xindex1 + 1, (int)c_histogramBucketCount - 1);
		float xindexfract = xindexf - std::floor(xindexf);

		float y1 = CDF[xindex1];
		float y2 = CDF[xindex2];

		csv[csvcolumnIndex + 1].values[index] = Lerp(y1, y2, xindexfract);
	}

	// Put the CDF into the CDF csv
	int cdfcsvcolumnIndex = (int)CDFcsv.size();
	CDFcsv.resize(cdfcsvcolumnIndex + 2);
	CDFcsv[cdfcsvcolumnIndex].label = std::string(label) + " CDF";
	CDFcsv[cdfcsvcolumnIndex].values = CDF;

	// Fit the CDF with a polynomial
	LeastSquaresPolynomialFit<3> fit;
	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float x = float(i) / float(CDF.size());
		float y = CDF[i];
		fit.AddPoint(x, y);
	}
	fit.CalculateCoefficients();

	// print out the coefficients
	printf("  y = ");
	bool first = true;
	for (size_t i = 0; i < fit.m_coefficients.size(); ++i)
	{
		printf("%s%f x^%i", first ? "" : " + ", fit.m_coefficients[fit.m_coefficients.size() - i - 1], (int)i);
		first = false;
	}
	printf("\n");

	// Put the values through the polynomial fit CDF (inverted, inverted CDF) to make them be a uniform distribution
	csv[csvcolumnIndex + 2].label = std::string(label) + "_ToUniformQuadraticFit";
	csv[csvcolumnIndex + 2].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[csvcolumnIndex].values[index];
		csv[csvcolumnIndex + 2].values[index] = fit.Evaluate(x);
	}

	// Put the fit CDF into the CDF csv
	CDFcsv[cdfcsvcolumnIndex + 1].label = std::string(label) + " Cubic";
	CDFcsv[cdfcsvcolumnIndex + 1].values.resize(CDF.size());
	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float x = float(i) / float(CDF.size() - 1);
		CDFcsv[cdfcsvcolumnIndex + 1].values[i] = fit.Evaluate(x);
	}
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

- try a piecewise fit least squares fit
 - could try a couple different degree fits and piecewise, and use whatever fits best.

NOTES:

- this has a more direct solve using gauss jordan, than inverting a matrix and multiplying
 - https://blog.demofox.org/2022/06/29/piecewise-least-squares-curve-fitting/
- could also mention the piecewise and weighted least squares fitting.

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
