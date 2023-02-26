#pragma once

#include <array>

template <size_t N>
class LeastSquaresPolynomialFit
{
public:

	void AddPoint(float x, float y)
	{
		double xpow = 1.0;
		for (int i = 0; i < m_ATA.size(); ++i)
		{
			m_ATA[i] += xpow;
			xpow *= x;
		}

		xpow = 1.0f;
		for (int i = 0; i < m_ATY.size(); ++i)
		{
			m_ATY[i] += xpow * y;
			xpow *= x;
		}
	}

	void CalculateCoefficients()
	{
		// make the full ATA matrix, as an augmeted matrix with the identity matrix on the right
		std::array<std::array<double, 2 * (N + 1)>, N + 1> ATA;
		for (int iy = 0; iy < N + 1; ++iy)
		{
			for (int ix = 0; ix < N + 1; ++ix)
			{
				ATA[iy][ix] = m_ATA[ix + iy];
				ATA[iy][ix + N + 1] = (ix == iy) ? 1.0 : 0.0;
			}
		}

		// invert the ATA matrix
		for (int targetColumn = 0; targetColumn < N + 1; ++targetColumn)
		{
			// Find the row with the biggest value in this column
			int bestRow = targetColumn;
			double bestRowValue = ATA[targetColumn][targetColumn];
			for (int row = targetColumn + 1; row < N + 1; ++row)
			{
				if (ATA[row][targetColumn] > bestRowValue)
				{
					bestRow = row;
					bestRowValue = ATA[row][targetColumn];
				}
			}

			// Divide the row by the value to make a 1 in the column
			for (int column = 0; column < 2 * (N + 1); ++column)
				ATA[bestRow][column] /= bestRowValue;

			// Subtract multiples of this row from other rows to make them have a 0 in this column
			for (int row = 0; row < N + 1; ++row)
			{
				if (row == bestRow)
					continue;

				double multiplier = ATA[row][targetColumn];

				for (int column = 0; column < 2 * (N + 1); ++column)
					ATA[row][column] -= ATA[bestRow][column] * multiplier;
			}

			// Swap this row into the targetColumn'th row
			if (bestRow != targetColumn)
			{
				for (int column = 0; column < 2 * (N + 1); ++column)
					std::swap(ATA[bestRow][column], ATA[targetColumn][column]);
			}
		}

		// Multiply ATA^-1 by ATY to get the coefficients
		for (int coefficientIndex = 0; coefficientIndex < N + 1; ++coefficientIndex)
		{
			m_coefficients[coefficientIndex] = 0.0f;
			for (int index = 0; index < N + 1; ++index)
				m_coefficients[coefficientIndex] += ATA[coefficientIndex][N + 1 + index] * m_ATY[index];
		}
	}

	float Evaluate(float x)
	{
		double ret = 0.0;

		double xpow = 1.0;
		for (int index = 0; index < N + 1; ++index)
		{
			ret += xpow * m_coefficients[index];
			xpow *= x;
		}

		return (float)ret;
	}

private:
	std::array<double, (N + 1) * 2 - 1> m_ATA = {};
	std::array<double, N + 1> m_ATY = {};

	std::array<double, N + 1> m_coefficients = {};
};
