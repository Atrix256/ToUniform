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

		// a point constraint between each piece
		// also a point constraint for f(0) = 0 and f(1) = 1
		static const int ATAMatrixConstraints = PIECES - 1 + 2;

		static const int ATAMatrixHeightNoConstraints = ATAPieceWidth * PIECES;
		static const int ATAMatrixHeight = ATAMatrixHeightNoConstraints + ATAMatrixConstraints;
		static const int ATAMtrixAugmentedWidth = ATAMatrixHeight + 1;
		static const int ATAMatrixCoefficients = ATAPieceWidth * PIECES;

		// make the ATA matrix, as an augmeted matrix with the ATY matrix on the right.
		std::array<std::array<double, ATAMtrixAugmentedWidth>, ATAMatrixHeight> ATA = {};
		for (int iy = 0; iy < ATAMatrixHeightNoConstraints; ++iy)
		{
			ATA[iy][ATAMtrixAugmentedWidth - 1] = m_ATY[iy / ATAPieceWidth][iy % ATAPieceWidth];

			int pieceIndex = iy / ATAPieceWidth;
			for (int ix = 0; ix < ATAPieceWidth; ++ix)
				ATA[iy][ix + pieceIndex * ATAPieceWidth] = m_ATA[pieceIndex][ix + iy % ATAPieceWidth];
		}

		// fill out the C0 constraints between pieces
		for (int constraint = 0; constraint < PIECES - 1; ++constraint)
		{
			int pieceIndex1 = constraint;
			int pieceIndex2 = constraint + 1;

			int row = ATAMatrixHeightNoConstraints + constraint;

			float x = float(pieceIndex2) / float(PIECES);

			double xpow = 1.0;
			for (int index = 0; index < ORDER + 1; ++index)
			{
				// C
				ATA[row][pieceIndex1 * (ORDER + 1) + index] = xpow;
				ATA[row][pieceIndex2 * (ORDER + 1) + index] = -xpow;

				// C^T
				ATA[pieceIndex1 * (ORDER + 1) + index][ATAPieceWidth * PIECES + constraint] = xpow;
				ATA[pieceIndex2 * (ORDER + 1) + index][ATAPieceWidth * PIECES + constraint] = -xpow;

				xpow *= x;
			}
		}

		// make constraints for f(0) = 0 and f(1) = 1
		for (int constraint = PIECES - 1; constraint < ATAMatrixConstraints; ++constraint)
		{
			int row = ATAMatrixHeightNoConstraints + constraint;

			float x = (constraint == PIECES - 1) ? 0.0f : 1.0f;
			float z = x;
			int pieceIndex = (constraint == PIECES - 1) ? 0 : PIECES - 1;

			// Z
			ATA[row][ATAMtrixAugmentedWidth - 1] = z;

			double xpow = 1.0;
			for (int index = 0; index < ORDER + 1; ++index)
			{
				// C
				ATA[row][pieceIndex * (ORDER + 1) + index] = xpow;

				// C^T
				ATA[pieceIndex * (ORDER + 1) + index][ATAPieceWidth * PIECES + constraint] = xpow;

				xpow *= x;
			}
		}

		// invert the ATA matrix
		for (int targetColumn = 0; targetColumn < ATAMtrixAugmentedWidth - 1; ++targetColumn)
		{
			// Find the row with the biggest value in this column
			int bestRow = targetColumn;
			double bestRowValue = ATA[targetColumn][targetColumn];
			for (int row = targetColumn + 1; row < ATAMtrixAugmentedWidth - 1; ++row)
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
			for (int row = 0; row < ATAMatrixHeight; ++row)
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
		for (int coefficientIndex = 0; coefficientIndex < ATAMatrixCoefficients; ++coefficientIndex)
			m_coefficients[coefficientIndex / ATAPieceWidth][coefficientIndex % ATAPieceWidth] = ATA[coefficientIndex][ATAMtrixAugmentedWidth - 1];
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
