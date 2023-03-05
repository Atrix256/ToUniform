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

std::string FindBestPolynomialFit(const std::vector<float>& CDF, CSV& csv, CSV& CDFcsv, int csvcolumnIndex, int cdfcsvcolumnIndex, const char* label)
{
	float lowestRMSE = FLT_MAX;
	std::string bestFormula;
	FindBestPolynomialFit_Order<1>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	FindBestPolynomialFit_Order<2>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	FindBestPolynomialFit_Order<3>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	//FindBestPolynomialFit_Order<4>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	//FindBestPolynomialFit_Order<5>(CDF, lowestRMSE, bestFormula, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
	printf("%s", bestFormula.c_str());
	return bestFormula;
}

void SequenceTest(CSV& csv, CSV& CDFcsv, int csvcolumnIndex, const char* label)
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

		// Make the error be split evenly on the positive and negative side
		for (size_t index = 1; index < CDFFull.size() - 1; ++index)
			CDFFull[index] = (CDFFull[index] + CDFFull[index + 1]) / 2.0f;
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

		// Make the error be split evenly on the positive and negative side
		for (size_t index = 1; index < CDFSmall.size() - 1; ++index)
			CDFSmall[index] = (CDFSmall[index] + CDFSmall[index + 1]) / 2.0f;
	}

	// Put the values through the full CDF (inverted, inverted CDF) to make them be a uniform distribution
	csv[csvcolumnIndex + 1].label = std::string(label) + "_ToUniform1024";
	csv[csvcolumnIndex + 1].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[csvcolumnIndex].values[index];

		float xindexf = std::min(x * float(CDFFull.size() - 1) - 0.5f, (float)(CDFFull.size() - 1));
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

		float xindexf = std::min(x * float(CDFSmall.size() - 1) - 0.5f, (float)(CDFSmall.size() - 1));
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
	CDFcsv[cdfcsvcolumnIndex].label = std::string(label) + " CDF 1024";
	CDFcsv[cdfcsvcolumnIndex].values = CDFFull;

	// expand the small cdf so that it looks right in the graphs
	{
		CDFcsv[cdfcsvcolumnIndex + 1].label = std::string(label) + " CDF 64";
		CDFcsv[cdfcsvcolumnIndex + 1].values.resize(CDFFull.size());
		for (size_t index = 0; index < CDFFull.size(); ++index)
		{
			float x = (float(index) + 0.5f) / float(CDFFull.size() - 1);
			int srcIndex = std::min(int(x * float(CDFSmall.size() - 1)), (int)CDFSmall.size() - 1);
			CDFcsv[cdfcsvcolumnIndex + 1].values[index] = CDFSmall[srcIndex];
		}
	}

	// Find the best piecewise polynomial fit we can for this CDF, and use that
	std::string bestFormula = FindBestPolynomialFit(CDFFull, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);

	// write to out.txt
	{
		FILE* file = nullptr;
		fopen_s(&file, "out.txt", "ab");

		fprintf(file, "\n==========================\n%s\n==========================\n\n", label);

		// write the polynomial
		fprintf(file, "%s\n", bestFormula.c_str());

		// write the small LUT
		fprintf(file, "float LUT[%i] = {\n", (int)CDFSmall.size());
		for (size_t i = 0; i < CDFSmall.size(); ++i)
		{
			if (i > 0 && (i % 10) == 0)
				fprintf(file, "    %f,  // %i\n", CDFSmall[i], (int)i);
			else
				fprintf(file, "    %f,\n", CDFSmall[i]);
		}
		fprintf(file, "};\n\n");

		fclose(file);
	}
}

void FinalBNTests(pcg32_random_t& rng, CSV& csv, CSV& CDFcsv)
{
	printf("\nFinal Blue Noise\n");

	// reserve space in the CSV for this data 
	// 0) The filtered white noise
	// 1) To uniform with 1024 table CDF
	// 2) To uniform with 64 table CDF
	// 3) To uniform with least squares
	int csvcolumnIndex1 = (int)csv.size();
	int csvcolumnIndex2 = csvcolumnIndex1 + 4;
	csv.resize(csv.size() + 4); // TODO: +8
	csv[csvcolumnIndex1].label = "Final BN LUT";
	//csv[csvcolumnIndex2].label = "Final BN Polynomial";

	// make white noise
	std::vector<float> whiteNoise(c_numberCount);
	for (float& f : whiteNoise)
		f = PCGRandomFloat01(rng);

	static const std::vector<float> xCoefficients = { 0.5f, -1.0f, 0.5f };
	static const std::vector<float> yCoefficients = { };

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

	// Convert to uniform using the LUT
	{
		csv[csvcolumnIndex1].values = filteredWhiteNoise;
		for (float& f : csv[csvcolumnIndex1].values)
		{
			// TODO: this
		}
	}

	// Convert to uniform using the polynomial
	{
		// TODO: temp
		csvcolumnIndex2 = csvcolumnIndex1;

		auto ToUniform = [](float x)
		{
			int bucket = std::min(int(x * 4.0f), 3);
			float ret = 0.0f;

			// the second half of the curve is the first curve, flipped on the x and y axis
			if (bucket & 2)
			{
				x = 1.0f - x;
				bucket ^= 1;
			}

			if ((bucket & 1) == 0)
				ret = 5.26063f * x * x * x + 0.0402881f * x * x + 0.000498876f * x + 0.0f;
			else
				ret = -5.21888f * x * x * x + 7.8326f * x * x - 1.93075f * x + 0.159535f;

			if (bucket & 2)
				ret = 1.0f - ret;

			return ret;
		};

		// Normalize to [0,1]
		float themin = filteredWhiteNoise[0];
		float themax = filteredWhiteNoise[0];
		for (float f : filteredWhiteNoise)
		{
			themin = std::min(themin, f);
			themax = std::max(themax, f);
		}
		for (float& f : filteredWhiteNoise)
			f = (f - themin) / (themax - themin);

		// convert to uniform
		csv[csvcolumnIndex2].values = filteredWhiteNoise;
		for (float& f : csv[csvcolumnIndex2].values)
			f = ToUniform(f);

		// do the normal tests
		SequenceTest(csv, CDFcsv, csvcolumnIndex2, "Final BN Polynomial");
	}

	// Do the rest of the testing
	//SequenceTest(csv, CDFcsv, csvcolumnIndex, label);

	/*
 Order 3 with 4 pieces. RMSE = 4.61208e-05
 x in [0, 0.25)
  y = 5.26063 x^3 + 0.0402881 x^2 + 0.000498876 x + 0
 x in [0.25, 0.5)
  y = -5.21888 x^3 + 7.8326 x^2 + -1.93075 x + 0.159535
 x in [0.5, 0.75)
  y = -5.27361 x^3 + 7.92327 x^2 + -1.98037 x + 0.168519
 x in [0.75, 1]
  y = 5.23633 x^3 + -15.7565 x^2 + 15.8038 x + -4.28361

float LUT[65] = {
	0.000000,
	0.000115,
	0.000417,
	0.001032,
	0.002075,
	0.003690,
	0.005997,
	0.009095,
	0.013108,
	0.018147,
	0.024355,  // 10
	0.031841,
	0.040706,
	0.051112,
	0.063178,
	0.076976,
	0.092619,
	0.110101,
	0.129387,
	0.150341,
	0.172740,  // 20
	0.196455,
	0.221397,
	0.247503,
	0.274657,
	0.302720,
	0.331587,
	0.361112,
	0.391223,
	0.421820,
	0.452615,  // 30
	0.483524,
	0.514501,
	0.545454,
	0.576250,
	0.606779,
	0.636905,
	0.666491,
	0.695483,
	0.723680,
	0.750922,  // 40
	0.777081,
	0.802167,
	0.826038,
	0.848501,
	0.869477,
	0.888797,
	0.906394,
	0.922151,
	0.936074,
	0.948246,  // 50
	0.958744,
	0.967731,
	0.975299,
	0.981571,
	0.986677,
	0.990725,
	0.993860,
	0.996189,
	0.997827,
	0.998915,  // 60
	0.999560,
	0.999873,
	0.999980,
	1.000000,
};
	*/
}

void VoidAndClusterTest(pcg32_random_t& rng, CSV& csv, CSV& CDFcsv)
{
	const char* label = "VoidAndCluster";
	printf("\n%s\n", label);

	// reserve space in the CSV for this data 
	// 0) The filtered white noise
	// 1) To uniform with 1024 table CDF
	// 2) To uniform with 64 table CDF
	// 3) To uniform with least squares
	int csvcolumnIndex = (int)csv.size();
	csv.resize(csvcolumnIndex + 4);
	csv[csvcolumnIndex].label = label;

	// read the data from disk
	std::vector<size_t> blueNoise;
	for (int i = 0; i < 100; ++i)
	{
		char fileName[256];
		sprintf_s(fileName, "bluenoise/bn100k_%i.bin", i);
		FILE* file = nullptr;
		fopen_s(&file, fileName, "rb");
		size_t length = 0;
		fread(&length, sizeof(size_t), 1, file);
		size_t oldLength = blueNoise.size();
		blueNoise.resize(oldLength + length);
		fread(&blueNoise[oldLength], sizeof(size_t), length, file);
		fclose(file);
	}

	size_t length = blueNoise.size();
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
	SequenceTest(csv, CDFcsv, csvcolumnIndex, label);
}

void IIRTest(const char* label, pcg32_random_t& rng, CSV& csv, CSV& CDFcsv, const std::vector<float>& xCoefficients, const std::vector<float>& yCoefficients)
{
	printf("\n%s\n", label);

	// reserve space in the CSV for this data 
	// 0) The filtered white noise
	// 1) To uniform with 1024 table CDF
	// 2) To uniform with 64 table CDF
	// 3) To uniform with least squares
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

	// reserve space in the CSV for this data 
	// 0) The filtered white noise
	// 1) To uniform with 1024 table CDF
	// 2) To uniform with 64 table CDF
	// 3) To uniform with least squares
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

	// empty out.txt
	{
		FILE* file = nullptr;
		fopen_s(&file, "out.txt", "wb");
		fclose(file);
	}

	// make noise and add to csv
	CSV csv, CDFcsv;
	/*
	FIRTest("Box3RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f });
	//FIRTest("Box3RedNoise2", rng, csv, CDFcsv, { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f });
	FIRTest("Box3BlueNoise", rng, csv, CDFcsv, { -1.0f, 1.0f, -1.0f });

	FIRTest("Box5RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f });
	FIRTest("Box5BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f, 1.0f, -1.0f, 1.0f });

	//FIRTest("Gauss10BlueNoiseNarrow", rng, csv, CDFcsv, { -0.0002f, -0.0060f, -0.0606f, -0.2417f, 1.0f - 0.3829f, -0.2417f, -0.0606f, -0.0060f, -0.0002f });

	FIRTest("Gauss10BlueNoise", rng, csv, CDFcsv, { 0.0002f, -0.0060f, 0.0606f, -0.2417f, 0.3829f, -0.2417f, 0.0606f, -0.0060f, 0.0002f });

	IIRTest("FIRHPF", rng, csv, CDFcsv, { 0.5f, -1.0f, 0.5f }, {});
	IIRTest("IIRHPF", rng, csv, CDFcsv, { 0.5f, -1.0f, 0.5f }, { 0.9f });

	VoidAndClusterTest(rng, csv, CDFcsv);
	*/

	FinalBNTests(rng, csv, CDFcsv);

	printf("\nWriting CSVs...\n");
	WriteCSV(csv, "out.csv");
	WriteCSV(CDFcsv, "cdf.csv");

	printf("\nRunning MakeHistograms.py\n");
	system("python MakeHistograms.py");

	return 0;
}

/*
TODO:

* make a copy/pastable class for the LUT and polynomial blue noise

! should have some simple code to make colored uniform noise by the end. need it for the next post!
 * FIRHPF put through the O3 C4 function
 ! the code should just have the first half of the curve, and mirror / reverse x and y to make the second half.
 * we can also make a 32 value LUT and use it the same way
 * put both of these methods through the tests above.

! need to fix or hide how uniform PDFs get mangled by the histogram etc code. Probably mark it as "don't test" or just fill in the data with what it should be after with a comment about it. and comment in blog post.

Alan: make sure folks understood how any stocahastic rendering can be made into an ML problem.
 * also initing with good starting values (noise).
 * stochastic rendering is anything that uses random values (scalars or vectors or whatever else).
 * noise is a "blind attempt" at making the result good.
 * ML at runtime makes that not blind.
 * algorithms: AO, dithering, rendering volumetric clouds. motion blur, DOF, material blending ....

NOTES:

* talk about all the point and derivative constraints
* talk about how you made the LUT (made [32] be 0.5) and the polynomial.

* Lut down to 64 entries wasn't bad.
 * probably would be better if the lut was non linear spaced points but that'd be hard to look up

* blue noise made by void and cluster, using a sigma of 1.0
 * blue noise "tiles well" so may be fine using a small amount of blue noise values (a few 100?) and re-using them.

* C1 continuity helped the "devil horns" smooth out a bit, but they are still there

FIR: 0.5, -2, 1
http://demofox.org/DSPFIR/FIR.html
y = 0.500 * x(n) - 1.000 * x(n-1) + 0.500 * x(n-2)

IIR: 0.5, -2, 1, -0.9, 0
http://demofox.org/DSPIIR/IIR.html
y = 0.500 * x(n) - 1.000 * x(n-1) + 0.500 * x(n-2) + 0.900 * y(n-1)

* This is all FIR filters.
 * IIR filters could be interesting to explore, for making a stream of colored noise, and make it uniform with CDFs / fit CDFs

* a LUT is basically a piecewise linear curve with C0.
 * could have it be piecewise quadratic with c1 as well.

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
