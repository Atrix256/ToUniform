#pragma once

class BlueNoiseStreamLUT
{
public:
	BlueNoiseStreamLUT(pcg32_random_t rng)
		: m_rng(rng)
	{
		m_lastValues[0] = RandomFloat01();
		m_lastValues[1] = RandomFloat01();
	}

	float Next()
	{
		// Filter uniform white noise to remove low frequencies and make it blue.
		// A side effect is the noise becomes non uniform.
		static const float xCoefficients[3] = {0.5f, -1.0f, 0.5f};

		float value = RandomFloat01();

		float y =
			value * xCoefficients[0] +
			m_lastValues[0] * xCoefficients[1] +
			m_lastValues[1] * xCoefficients[2];

		m_lastValues[1] = m_lastValues[0];
		m_lastValues[0] = value;

		// the noise is also [-1,1] now, normalize to [0,1]
		float x = y * 0.5f + 0.5f;

		// Make the noise uniform again by putting it through a LUT of the CDF
		static const float LUT[] =
		{
			0.000008f,
			0.000126f,
			0.000478f,
			0.001155f,
			0.002251f,
			0.004030f,
			0.006441f,
			0.009596f,
			0.013695f,
			0.018927f,
			0.025355f,  // 10
			0.033032f,
			0.042107f,
			0.052496f,
			0.064787f,
			0.078712f,
			0.094684f,
			0.112292f,
			0.131824f,
			0.152818f,
			0.175335f,  // 20
			0.199049f,
			0.224476f,
			0.250783f,
			0.277852f,
			0.305973f,
			0.334680f,
			0.363986f,
			0.393945f,
			0.424846f,
			0.455969f,  // 30
			0.486683f,
			0.517827f,
			0.548521f,
			0.579180f,
			0.609718f,
			0.639919f,
			0.669180f,
			0.698089f,
			0.726230f,
			0.753389f,  // 40
			0.779630f,
			0.804358f,
			0.828272f,
			0.850459f,
			0.871144f,
			0.890308f,
			0.907814f,
			0.923487f,
			0.937238f,
			0.949365f,  // 50
			0.959586f,
			0.968409f,
			0.975893f,
			0.982001f,
			0.987002f,
			0.990933f,
			0.994027f,
			0.996302f,
			0.997890f,
			0.998971f,  // 60
			0.999558f,
			0.999887f,
			0.999991f
		};
		static const size_t lutSize = _countof(LUT);

		float xindexf = std::min(x * float(lutSize - 1), (float)(lutSize - 1));
		int xindex1 = int(xindexf);
		int xindex2 = std::min(xindex1 + 1, (int)lutSize - 1);
		float xindexfract = xindexf - std::floor(xindexf);

		float y1 = LUT[xindex1];
		float y2 = LUT[xindex2];

		if (xindex1 == 0 && xindexfract == 0.0f)
			y1 = y2 = 0.0f;
		else if (xindex1 == lutSize - 1)
			y1 = y2 = 1.0f;

		return Lerp(y1, y2, xindexfract);
	}

private:
	float RandomFloat01()
	{
		// return a uniform white noise random float between 0 and 1.
		// Can use whatever RNG you want, such as std::mt19937.
		return ldexpf((float)pcg32_random_r(&m_rng), -32);
	}

	pcg32_random_t m_rng;
	float m_lastValues[2] = {};
};

class BlueNoiseStreamPolynomial
{
public:
	BlueNoiseStreamPolynomial(pcg32_random_t rng)
		: m_rng(rng)
	{
		m_lastValues[0] = RandomFloat01();
		m_lastValues[1] = RandomFloat01();
	}
	
	float Next()
	{
		// Filter uniform white noise to remove low frequencies and make it blue.
		// A side effect is the noise becomes non uniform.
		static const float xCoefficients[3] = {0.5f, -1.0f, 0.5f};

		float value = RandomFloat01();

		float y =
			value * xCoefficients[0] +
			m_lastValues[0] * xCoefficients[1] +
			m_lastValues[1] * xCoefficients[2];

		m_lastValues[1] = m_lastValues[0];
		m_lastValues[0] = value;

		// the noise is also [-1,1] now, normalize to [0,1]
		float x = y * 0.5f + 0.5f;

		// Make the noise uniform again by putting it through a piecewise cubic polynomial approximation of the CDF
		// Switched to Horner's method polynomials, and a polynomial array to avoid branching, per Marc Reynolds. Thanks!
		float polynomialCoefficients[16] = {
			5.25964f, 0.039474f, 0.000708779f, 0.0f,
			-5.20987f, 7.82905f, -1.93105f, 0.159677f,
			-5.22644f, 7.8272f, -1.91677f, 0.15507f,
			5.23882f, -15.761f, 15.8054f, -4.28323f
		};
		int first = std::min(int(x * 4.0f), 3) * 4;
		return polynomialCoefficients[first + 3] + x * (polynomialCoefficients[first + 2] + x * (polynomialCoefficients[first + 1] + x * polynomialCoefficients[first + 0]));
	}

private:
	float RandomFloat01()
	{
		// return a uniform white noise random float between 0 and 1.
		// Can use whatever RNG you want, such as std::mt19937.
		return ldexpf((float)pcg32_random_r(&m_rng), -32);
	}

	pcg32_random_t m_rng;
	float m_lastValues[2] = {};
};

class RedNoiseStreamPolynomial
{
public:
	RedNoiseStreamPolynomial(pcg32_random_t rng)
		: m_rng(rng)
	{
		m_lastValues[0] = RandomFloat01();
		m_lastValues[1] = RandomFloat01();
	}

	float Next()
	{
		// Filter uniform white noise to remove high frequencies and make it red.
		// A side effect is the noise becomes non uniform.
		static const float xCoefficients[3] = { 0.25f, 0.5f, 0.25f };

		float value = RandomFloat01();

		float y =
			value * xCoefficients[0] +
			m_lastValues[0] * xCoefficients[1] +
			m_lastValues[1] * xCoefficients[2];

		m_lastValues[1] = m_lastValues[0];
		m_lastValues[0] = value;

		float x = y;

		// Make the noise uniform again by putting it through a piecewise cubic polynomial approximation of the CDF
		// Switched to Horner's method polynomials, and a polynomial array to avoid branching, per Marc Reynolds. Thanks!
		float polynomialCoefficients[16] = {
			5.25964f, 0.039474f, 0.000708779f, 0.0f,
			-5.20987f, 7.82905f, -1.93105f, 0.159677f,
			-5.22644f, 7.8272f, -1.91677f, 0.15507f,
			5.23882f, -15.761f, 15.8054f, -4.28323f
		};
		int first = std::min(int(x * 4.0f), 3) * 4;
		return polynomialCoefficients[first + 3] + x * (polynomialCoefficients[first + 2] + x * (polynomialCoefficients[first + 1] + x * polynomialCoefficients[first + 0]));
	}

private:
	float RandomFloat01()
	{
		// return a uniform white noise random float between 0 and 1.
		// Can use whatever RNG you want, such as std::mt19937.
		return ldexpf((float)pcg32_random_r(&m_rng), -32);
	}

	pcg32_random_t m_rng;
	float m_lastValues[2] = {};
};