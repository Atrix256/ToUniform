#pragma once

#include <array>

template <size_t ORDER, size_t PIECES>
class LeastSquaresPolynomialFit
{
public:

	void AddPoint(float x, float y)
	{
		int bucket = std::min(int(x * float(PIECES)), (int)PIECES - 1);

		double xpow = 1.0;
		for (int i = 0; i < m_ATA[bucket].size(); ++i)
		{
			m_ATA[bucket][i] += xpow;
			xpow *= x;
		}

		xpow = 1.0;
		for (int i = 0; i < m_ATY[bucket].size(); ++i)
		{
			m_ATY[bucket][i] += xpow * y;
			xpow *= x;
		}
	}

	void CalculateCoefficients()
	{
		static const int ATAPieceWidth = ORDER + 1;
		static const int ATAMatrixWidthHeight = ATAPieceWidth * PIECES;
		static const int ATAMtrixAugmentedWidth = ATAMatrixWidthHeight + 1;

		// make the full ATA matrix, as an augmeted matrix with the ATY matrix on the right
		std::array<std::array<double, ATAMtrixAugmentedWidth>, ATAMatrixWidthHeight> ATA = {};
		for (int iy = 0; iy < ATAMatrixWidthHeight; ++iy)
		{
			ATA[iy][ATAMatrixWidthHeight] = m_ATY[iy / ATAPieceWidth][iy % ATAPieceWidth];

			int pieceIndex = iy / ATAPieceWidth;
			for (int ix = 0; ix < ATAPieceWidth; ++ix)
				ATA[iy][ix + pieceIndex * ATAPieceWidth] = m_ATA[pieceIndex][ix + iy % ATAPieceWidth];
		}

		// invert the ATA matrix
		for (int targetColumn = 0; targetColumn < ATAMatrixWidthHeight; ++targetColumn)
		{
			// Find the row with the biggest value in this column
			int bestRow = targetColumn;
			double bestRowValue = ATA[targetColumn][targetColumn];
			for (int row = targetColumn + 1; row < ATAMatrixWidthHeight; ++row)
			{
				if (std::abs(ATA[row][targetColumn]) > std::abs(bestRowValue))
				{
					bestRow = row;
					bestRowValue = ATA[row][targetColumn];
				}
			}

			// Divide the row by the value to make a 1 in the column
			for (int column = 0; column < ATAMtrixAugmentedWidth; ++column)
				ATA[bestRow][column] /= bestRowValue;

			// Subtract multiples of this row from other rows to make them have a 0 in this column
			for (int row = 0; row < ATAMatrixWidthHeight; ++row)
			{
				if (row == bestRow)
					continue;

				double multiplier = ATA[row][targetColumn];

				for (int column = 0; column < ATAMtrixAugmentedWidth; ++column)
					ATA[row][column] -= ATA[bestRow][column] * multiplier;
			}

			// Swap this row into the targetColumn'th row
			if (bestRow != targetColumn)
			{
				for (int column = 0; column < ATAMtrixAugmentedWidth; ++column)
					std::swap(ATA[bestRow][column], ATA[targetColumn][column]);
			}
		}

		// The coefficients are on the right side
		for (int coefficientIndex = 0; coefficientIndex < ATAMatrixWidthHeight; ++coefficientIndex)
			m_coefficients[coefficientIndex / ATAPieceWidth][coefficientIndex % ATAPieceWidth] = ATA[coefficientIndex][ATAMatrixWidthHeight];
	}

	float Evaluate(float x)
	{
		int bucket = std::min(int(x * float(PIECES)), (int)PIECES - 1);

		double ret = 0.0;

		double xpow = 1.0;
		for (int index = 0; index < ORDER + 1; ++index)
		{
			ret += xpow * m_coefficients[bucket][index];
			xpow *= x;
		}

		return (float)ret;
	}

public:
	std::array<std::array<double, ORDER + 1>, PIECES> m_coefficients = {};

private:
	std::array<std::array<double, (ORDER + 1) * 2 - 1>, PIECES> m_ATA = {};
	std::array<std::array<double, ORDER + 1>, PIECES> m_ATY = {};
};
