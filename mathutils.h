#pragma once

#include <vector>

inline float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

inline std::vector<float> Convolve(const std::vector<float>& A, const std::vector<float>& B)
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
