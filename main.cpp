#include <random>
#include <vector>
#include "pcg/pcg_basic.h"
#include <string>
#include "csv.h"
#include "mathutils.h"
#include "leastsquaresfit.h"
#include <sstream>

#define DETERMINISTIC() false

// The size of the list of random numbers output
static const size_t c_numberCount = 10000000;

// Bucket count of histogram that makes the PDF and CDF
static const size_t c_histogramBucketCountFull = 1024;
static const size_t c_histogramBucketCountSmall = 64;

float PCGRandomFloat01(pcg32_random_t& rng)
{
	return ldexpf((float)pcg32_random_r(&rng), -32);
}

template <size_t ORDER, size_t PIECES>
void FindBestPolynomialFit_Order_Pieces(const std::vector<float>& CDF, float& lowestRMSE, std::string& bestFormula, CSV& csv, CSV& CDFcsv, int csvcolumnIndex, int cdfcsvcolumnIndex, const char* label)
{
	// fit a piecewise polynomial to the CDF
	LeastSquaresPolynomialFit<ORDER, PIECES> fit;
	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float x = float(i) / float(CDF.size() - 1.0f);
		float y = CDF[i];
		fit.AddPoint(x, y);
	}

	fit.CalculateCoefficients();

	// Get the RMSE of this fit
	float RMSE = 0.0f;
	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float percent = float(i) / float(CDF.size() - 1);
		float error = CDF[i] - fit.Evaluate(percent);
		RMSE = Lerp(RMSE, error * error, 1.0f / float(i + 1));
	}
	RMSE = std::sqrt(RMSE);

	// if the RMSE is higher, don't take this fit
	if (RMSE >= lowestRMSE)
		return;

	// write the function so the best one can be printed out later
	std::stringstream formula;
	formula << " Order " << ORDER << " with " << PIECES << " pieces. RMSE = " << RMSE << "\n";
	for (size_t pieceIndex = 0; pieceIndex < PIECES; ++pieceIndex)
	{
		if (PIECES > 1)
		{
			float xmin = float(pieceIndex) / float(PIECES);
			float xmax = float(pieceIndex + 1) / float(PIECES);

			if (pieceIndex + 1 < PIECES)
				formula << " x in [" << xmin << ", " << xmax << ")\n";
			else
				formula << " x in [" << xmin << ", " << xmax << "]\n";
		}

		formula << "  y = ";
		bool first = true;
		for (size_t i = 0; i < fit.m_coefficients[pieceIndex].size(); ++i)
		{
			if (!first)
				formula << " + ";

			int xpower = int(fit.m_coefficients[pieceIndex].size() - i - 1);

			if (xpower == 0)
				formula << fit.m_coefficients[pieceIndex][xpower];
			else if (xpower == 1)
				formula << fit.m_coefficients[pieceIndex][xpower] << " x";
			else
				formula << fit.m_coefficients[pieceIndex][xpower] << " x^" << xpower;

			first = false;
		}
		formula << "\n";
	}
	bestFormula = formula.str();

	// Set the label
	char buffer[1024];
	sprintf_s(buffer, "%s_ToUniformFit_O%i_C%i", label, (int)ORDER, (int)PIECES);
	csv[csvcolumnIndex + 3].label = buffer;

	// Put the values through the polynomial fit CDF (inverted, inverted CDF) to make them be a uniform distribution
	csv[csvcolumnIndex + 3].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[csvcolumnIndex].values[index];
		csv[csvcolumnIndex + 3].values[index] = fit.Evaluate(x);
	}

	// Put the fit CDF into the CDF csv
	CDFcsv[cdfcsvcolumnIndex + 2].label = buffer;
	CDFcsv[cdfcsvcolumnIndex + 2].values.resize(CDF.size());
	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float x = float(i) / float(CDF.size() - 1);
		CDFcsv[cdfcsvcolumnIndex + 2].values[i] = fit.Evaluate(x);
	}
}

template <size_t ORDER>
void FindBestPolynomialFit_Order(const std::vector<float>& CDF, float& lowestRMSE, std::string& bestFormula, CSV& csv, CSV& CDFcsv, int csvcolumnIndex, int cdfcsvcolumnIndex, const char* label)
{
	FindBestPolynomialFit_Order_Pieces<ORDER, 1>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	FindBestPolynomialFit_Order_Pieces<ORDER, 2>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	FindBestPolynomialFit_Order_Pieces<ORDER, 3>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	FindBestPolynomialFit_Order_Pieces<ORDER, 4>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	//FindBestPolynomialFit_Order_Pieces<ORDER, 5>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
}

void FindBestPolynomialFit(const std::vector<float>& CDF, CSV& csv, CSV& CDFcsv, int csvcolumnIndex, int cdfcsvcolumnIndex, const char* label)
{
	float lowestRMSE = FLT_MAX;
	std::string bestFormula;
	FindBestPolynomialFit_Order<1>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	FindBestPolynomialFit_Order<2>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	FindBestPolynomialFit_Order<3>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	//FindBestPolynomialFit_Order<4>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	//FindBestPolynomialFit_Order<5>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	printf("%s", bestFormula.c_str());
}

void SequenceTest(CSV& csv, CSV& CDFcsv, int csvcolumnIndex, const char* label, bool isUniform = false)
{
	// Normalize it to [0,1] and put it into the csv
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
	std::vector<int> countsFull(c_histogramBucketCountFull, 0);
	std::vector<int> countsSmall(c_histogramBucketCountSmall, 0);
	for (float f : csv[csvcolumnIndex].values)
	{
		{
			int bucket = std::min(int(f * float(c_histogramBucketCountFull)), (int)c_histogramBucketCountFull - 1);
			countsFull[bucket]++;
		}
		{
			int bucket = std::min(int(f * float(c_histogramBucketCountSmall)), (int)c_histogramBucketCountSmall - 1);
			countsSmall[bucket]++;
		}
	}

	std::vector<float> CDFFull(c_histogramBucketCountFull);
	std::vector<float> CDFSmall(c_histogramBucketCountSmall);

	if (isUniform)
	{
		for (size_t index = 0; index < CDFFull.size(); ++index)
			CDFFull[index] = (float(index) + 0.5f) / float(CDFFull.size() - 1);

		for (size_t index = 0; index < CDFSmall.size(); ++index)
			CDFSmall[index] = (float(index) + 0.5f) / float(CDFSmall.size() - 1);
	}
	else
	{
		// Make Full CDF
		{
			std::vector<float> PDF(c_histogramBucketCountFull);
			float lastCDFValue = 0.0f;
			for (int index = 0; index < c_histogramBucketCountFull; ++index)
			{
				PDF[index] = float(countsFull[index]) / float(c_numberCount);
				CDFFull[index] = lastCDFValue + PDF[index];
				lastCDFValue = CDFFull[index];
			}
			CDFFull.insert(CDFFull.begin(), 0.0f);
			CDFFull[CDFFull.size() - 1] = 1.0f;
		}

		// Make Small CDF
		{
			std::vector<float> PDF(c_histogramBucketCountSmall);
			float lastCDFValue = 0.0f;
			for (int index = 0; index < c_histogramBucketCountSmall; ++index)
			{
				PDF[index] = float(countsSmall[index]) / float(c_numberCount);
				CDFSmall[index] = lastCDFValue + PDF[index];
				lastCDFValue = CDFSmall[index];
			}
			CDFSmall.insert(CDFSmall.begin(), 0.0f);
			CDFSmall[CDFSmall.size() - 1] = 1.0f;
		}
	}

	// Put the values through the full CDF (inverted, inverted CDF) to make them be a uniform distribution
	csv[csvcolumnIndex + 1].label = std::string(label) + "_ToUniform1024";
	csv[csvcolumnIndex + 1].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[csvcolumnIndex].values[index];

		float xindexf = std::min(x * float(CDFFull.size() - 1), (float)(CDFFull.size() - 1));
		int xindex1 = int(xindexf);
		int xindex2 = std::min(xindex1 + 1, (int)CDFFull.size() - 1);
		float xindexfract = xindexf - std::floor(xindexf);

		float y1 = CDFFull[xindex1];
		float y2 = CDFFull[xindex2];

		csv[csvcolumnIndex + 1].values[index] = Lerp(y1, y2, xindexfract);
	}

	// Put the values through the small CDF (inverted, inverted CDF) to make them be a uniform distribution
	csv[csvcolumnIndex + 2].label = std::string(label) + "_ToUniform64";
	csv[csvcolumnIndex + 2].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[csvcolumnIndex].values[index];

		float xindexf = std::min(x * float(CDFSmall.size() - 1), (float)(CDFSmall.size() - 1));
		int xindex1 = int(xindexf);
		int xindex2 = std::min(xindex1 + 1, (int)CDFSmall.size() - 1);
		float xindexfract = xindexf - std::floor(xindexf);

		float y1 = CDFSmall[xindex1];
		float y2 = CDFSmall[xindex2];

		csv[csvcolumnIndex + 2].values[index] = Lerp(y1, y2, xindexfract);
	}

	// Put the CDF into the CDF csv
	int cdfcsvcolumnIndex = (int)CDFcsv.size();
	CDFcsv.resize(cdfcsvcolumnIndex + 3);
	CDFcsv[cdfcsvcolumnIndex].label = std::string(label) + " CDF";
	CDFcsv[cdfcsvcolumnIndex].values = CDFFull;

	// expand the small cdf so that it looks right in the graphs
	{
		CDFcsv[cdfcsvcolumnIndex + 1].label = std::string(label) + " CDFSmall";
		CDFcsv[cdfcsvcolumnIndex + 1].values.resize(CDFFull.size());
		for (size_t index = 0; index < CDFFull.size(); ++index)
		{
			float x = (float(index) + 0.5f) / float(CDFFull.size() - 1);
			int srcIndex = std::min(int(x * float(CDFSmall.size() - 1)), (int)CDFSmall.size() - 1);
			CDFcsv[cdfcsvcolumnIndex + 1].values[index] = CDFSmall[srcIndex];
		}
	}

	// Find the best piecewise polynomial fit we can for this CDF, and use that
	FindBestPolynomialFit(CDFFull, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
}

void VoidAndClusterTest(pcg32_random_t& rng, CSV& csv, CSV& CDFcsv)
{
	const char* label = "VoidAndCluster";
	printf("\n%s\n", label);

	// reserve space in the CSV for this data.
	// 0) The noise from disk
	// 1) To uniform with CDF
	// 2) To uniform with least squares
	int csvcolumnIndex = (int)csv.size();
	csv.resize(csvcolumnIndex + 4);
	csv[csvcolumnIndex].label = label;

	// read the data from disk
	FILE* file = nullptr;
	fopen_s(&file, "bluenoise_100k_u64.bin", "rb");
	size_t length = 0;
	fread(&length, sizeof(size_t), 1, file);
	std::vector<size_t> blueNoise(length);
	fread(blueNoise.data(), sizeof(size_t), length, file);
	fclose(file);

	if (length > c_numberCount)
	{
		blueNoise.resize(c_numberCount);
		length = c_numberCount;
	}

	// convert to float
	csv[csvcolumnIndex].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
		csv[csvcolumnIndex].values[index] = float(blueNoise[index % length]) / float(length - 1);

	// Do the rest of the testing
	SequenceTest(csv, CDFcsv, csvcolumnIndex, label, true);
}

void IIRTest(const char* label, pcg32_random_t& rng, CSV& csv, CSV& CDFcsv, const std::vector<float>& xCoefficients, const std::vector<float>& yCoefficients)
{
	printf("\n%s\n", label);

	// reserve space in the CSV for this data.
	// 0) The filtered white noise
	// 1) To uniform with CDF
	// 2) To uniform with least squares
	int csvcolumnIndex = (int)csv.size();
	csv.resize(csvcolumnIndex + 4);
	csv[csvcolumnIndex].label = label;

	// make white noise
	std::vector<float> whiteNoise(c_numberCount);
	for (float& f : whiteNoise)
		f = PCGRandomFloat01(rng);

	// filter the white noise
	std::vector<float> filteredWhiteNoise(c_numberCount, 0.0f);
	for (int index = 0; index < c_numberCount; ++index)
	{
		float& out = filteredWhiteNoise[index];

		// FIR filtering
		for (size_t xIndex = 0; xIndex < xCoefficients.size(); ++xIndex)
		{
			if (xIndex > index)
				break;

			out += xCoefficients[xIndex] * whiteNoise[index - xIndex];
		}

		// IIR filtering feedback
		for (size_t yIndex = 0; yIndex < yCoefficients.size(); ++yIndex)
		{
			if (yIndex >= index)
				break;

			out += yCoefficients[yIndex] * filteredWhiteNoise[index - yIndex - 1];
		}
	}
	csv[csvcolumnIndex].values = filteredWhiteNoise;

	// Do the rest of the testing
	SequenceTest(csv, CDFcsv, csvcolumnIndex, label);
}

void FIRTest(const char* label, pcg32_random_t& rng, CSV& csv, CSV& CDFcsv, const std::vector<float>& kernel)
{
	printf("\n%s\n", label);

	// reserve space in the CSV for this data.
	// 0) The filtered white noise
	// 1) To uniform with CDF
	// 2) To uniform with least squares
	int csvcolumnIndex = (int)csv.size();
	csv.resize(csvcolumnIndex + 4);
	csv[csvcolumnIndex].label = label;

	// make white noise
	std::vector<float> whiteNoise(c_numberCount);
	for (float& f : whiteNoise)
		f = PCGRandomFloat01(rng);

	// Convolve and truncate the noise
	csv[csvcolumnIndex].values = Convolve(whiteNoise, kernel);
	csv[csvcolumnIndex].values.resize(c_numberCount);

	// Do the rest of the testing
	SequenceTest(csv, CDFcsv, csvcolumnIndex, label);
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
	FIRTest("Box3RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f });
	FIRTest("Box3RedNoise2", rng, csv, CDFcsv, { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f });
	FIRTest("Box3BlueNoise", rng, csv, CDFcsv, { -1.0f, 1.0f, -1.0f });

	FIRTest("Box5RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f });
	FIRTest("Box5BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f, 1.0f, -1.0f, 1.0f });

	//FIRTest("Gauss10BlueNoiseNarrow", rng, csv, CDFcsv, { -0.0002f, -0.0060f, -0.0606f, -0.2417f, 1.0f - 0.3829f, -0.2417f, -0.0606f, -0.0060f, -0.0002f });

	FIRTest("Gauss10BlueNoise", rng, csv, CDFcsv, { 0.0002f, -0.0060f, 0.0606f, -0.2417f, 0.3829f, -0.2417f, 0.0606f, -0.0060f, 0.0002f });

	IIRTest("FIRHPF", rng, csv, CDFcsv, { 0.5f, -1.0f, 0.5f }, {});
	IIRTest("IIRHPF", rng, csv, CDFcsv, { 0.5f, -1.0f, 0.5f }, { 0.9f });

	VoidAndClusterTest(rng, csv, CDFcsv);

	printf("\nWriting CSVs...\n");
	WriteCSV(csv, "out.csv");
	WriteCSV(CDFcsv, "cdf.csv");

	return 0;
}

/*
TODO:

* try c1 continuity? to help the spikes?
* maybe need c1 continuity too?
 * this basically could become a quadratic lut

* the lut should have the steps centered on the graph, not be on one side of it (add a half or something)
 * the uniform distribution has this, but not the other code path. it should!

* maybe see how low of a LUT size you can get away with?
 * down to 64 entries wasn't bad.
 * probably would be better if the lut was non linear spaced points but that'd be hard to look up
 ! for CDF graph, have a "high resolution LUT" version, and then a smaller bucket count verison too

! should have some simple code to make colored uniform noise by the end. need it for the next post!

NOTES:

* blue noise made by void and cluster, using a sigma of 1.0
 * blue noise "tiles well" so may be fine using a small amount of blue noise values (a few 100?) and re-using them.

FIR: 0.5, -2, 1
http://demofox.org/DSPFIR/FIR.html
y = 0.500 * x(n) - 1.000 * x(n-1) + 0.500 * x(n-2)

IIR: 0.5, -2, 1, -0.9, 0
http://demofox.org/DSPIIR/IIR.html
y = 0.500 * x(n) - 1.000 * x(n-1) + 0.500 * x(n-2) + 0.900 * y(n-1)

* This is all FIR filters.
 * IIR filters could be interesting to explore, for making a stream of colored noise, and make it uniform with CDFs / fit CDFs

* a LUT is basically a piecewise linear curve with C0.

? when making an LPF vs HPF, how do you control whether it's concave or convex?
 * well, this thread has great info: https://mastodon.gamedev.place/@demofox/109935390123342971

- this has a more direct solve using gauss jordan, than inverting a matrix and multiplying
 - https://blog.demofox.org/2022/06/29/piecewise-least-squares-curve-fitting/
- could also mention the piecewise and weighted least squares fitting.

- the best fit to the CDF seems to have the most segments and the highest order.
 - i limited it to 3 for both, but you could let it go higher in either
 - same shape histogram for same number of dice, so don't need to do this
 ? is that true for weighted dice though? no...

 * FIR calculator: http://demofox.org/DSPFIR/FIR.html
 * IIR calculator: http://demofox.org/DSPIIR/IIR.html

 ! Quick (?) blog post - better sharpening kernel by alternating pos and negative in low pass kernel.

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
