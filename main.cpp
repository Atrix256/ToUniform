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
static const size_t c_histogramBucketCount = 1024;

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
		float x = float(i) / float(CDF.size());
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
	csv[csvcolumnIndex + 2].label = buffer;

	// Put the values through the polynomial fit CDF (inverted, inverted CDF) to make them be a uniform distribution
	csv[csvcolumnIndex + 2].values.resize(c_numberCount);
	for (size_t index = 0; index < c_numberCount; ++index)
	{
		float x = csv[csvcolumnIndex].values[index];
		csv[csvcolumnIndex + 2].values[index] = fit.Evaluate(x);
	}

	// Put the fit CDF into the CDF csv
	CDFcsv[cdfcsvcolumnIndex + 1].label = buffer;
	CDFcsv[cdfcsvcolumnIndex + 1].values.resize(CDF.size());
	for (size_t i = 0; i < CDF.size(); ++i)
	{
		float x = float(i) / float(CDF.size() - 1);
		CDFcsv[cdfcsvcolumnIndex + 1].values[i] = fit.Evaluate(x);
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

void KernelTest(const char* label, pcg32_random_t& rng, CSV& csv, CSV& CDFcsv, const std::vector<float>& kernel)
{
	printf("\n%s\n", label);

	// reserve space in the CSV for this data.
	// 0) The filtered white noise
	// 1) To uniform with CDF
	// 2) To uniform with least squares
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

	// Find the best piecewise polynomial fit we can for this CDF, and use that
	FindBestPolynomialFit(CDF, csv, CDFcsv, csvcolumnIndex, cdfcsvcolumnIndex, label);
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
	KernelTest("Box4RedNoise", rng, csv, CDFcsv, { 1.0f, 1.0f, 1.0f, 1.0f });
	KernelTest("Box2BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f });
	KernelTest("Box4BlueNoise", rng, csv, CDFcsv, { 1.0f, -1.0f, 1.0f, -1.0f });

	KernelTest("Gauss10RedNoise", rng, csv, CDFcsv, { 0.0002f, 0.0060f, 0.0606f, 0.2417f, 0.3829f, 0.2417f, 0.0606f, 0.0060f, 0.0002f });

	KernelTest("Gauss10BlueNoise", rng, csv, CDFcsv, { -0.0002f, -0.0060f, -0.0606f, -0.2417f, 2.0f - 0.3829f, -0.2417f, -0.0606f, -0.0060f, -0.0002f });

	// TODO: more types of noise, and gaussian filtered.

	printf("\nWriting CSVs...\n");
	WriteCSV(csv, "out.csv");
	WriteCSV(CDFcsv, "cdf.csv");

	return 0;
}

/*
TODO:

* maybe see how low of a LUT size you can get away with?
* maybe need c1 continuity too? could also maybe have 0 derivatives on the sides?

* more noise types, including gaussian

! can we put a title for a figure above all the sub figures?

! histogram using fit curves doesn't seem right

! could maybe make a graph that showed the error between the two CDFs instead of showing them both next to each other.

NOTES:

- this has a more direct solve using gauss jordan, than inverting a matrix and multiplying
 - https://blog.demofox.org/2022/06/29/piecewise-least-squares-curve-fitting/
- could also mention the piecewise and weighted least squares fitting.

- the best fit to the CDF seems to have the most segments and the highest order.
 - i limited it to 3 for both, but you could let it go higher in either
 - same shape histogram for same number of dice, so don't need to do this
 ? is that true for weighted dice though? no...

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
