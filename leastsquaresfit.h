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

		// C0 constraints:
		//  2: f(0) = 0 and f(1) = 1
		//  PIECES - 1 : between each piece.
		// C1 constraints:
		//  PIECES - 1 : if ORDER > 1
		static const int PositionalConstraints = 2 + PIECES - 1;
		static const int DerivativeConstraints = (ORDER > 1) ? PIECES - 1 : 0;
		static const int ATAMatrixConstraints = PositionalConstraints + DerivativeConstraints;

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

		// The C0 constraints for f(0) = 0 and f(1) = 1
		for (int constraint = 0; constraint < 2; ++constraint)
		{
			int row = ATAMatrixHeightNoConstraints + constraint;

			float x = (constraint == 0) ? 0.0f : 1.0f;
			float z = (constraint == 0) ? 0.0f : 1.0f;
			int pieceIndex = (constraint == 0) ? 0 : PIECES - 1;

			// Z
			ATA[row][ATAMtrixAugmentedWidth - 1] = z;

			double xpow = 1.0;
			for (int index = 0; index < ATAPieceWidth; ++index)
			{
				// C
				ATA[row][pieceIndex * ATAPieceWidth + index] = xpow;

				// C^T
				ATA[pieceIndex * ATAPieceWidth + index][row] = xpow;

				xpow *= x;
			}
		}

		// The C0 constraints between pieces
		for (int constraint = 0; constraint < PIECES - 1; ++constraint)
		{
			int pieceIndex1 = constraint;
			int pieceIndex2 = constraint + 1;

			int row = ATAMatrixHeightNoConstraints + constraint + 2;

			float x = float(pieceIndex2) / float(PIECES);

			double xpow = 1.0;
			for (int index = 0; index < ATAPieceWidth; ++index)
			{
				// C
				ATA[row][pieceIndex1 * ATAPieceWidth + index] = xpow;
				ATA[row][pieceIndex2 * ATAPieceWidth + index] = -xpow;

				// C^T
				ATA[pieceIndex1 * ATAPieceWidth + index][row] = xpow;
				ATA[pieceIndex2 * ATAPieceWidth + index][row] = -xpow;

				xpow *= x;
			}
		}

		// The C1 constraints between pieces
		if (DerivativeConstraints > 0)
		{
			for (int constraint = 0; constraint < PIECES - 1; ++constraint)
			{
				int pieceIndex1 = constraint;
				int pieceIndex2 = constraint + 1;

				int row = ATAMatrixHeightNoConstraints + constraint + 2 + PIECES - 1;

				float x = float(pieceIndex2) / float(PIECES);

				double xpow = 1.0;
				for (int index = 0; index < ATAPieceWidth; ++index)
				{
					double value = double(index + 1) * xpow;

					// C
					ATA[row][pieceIndex1 * ATAPieceWidth + index] = value;
					ATA[row][pieceIndex2 * ATAPieceWidth + index] = -value;

					// C^T
					ATA[pieceIndex1 * ATAPieceWidth + index][row] = value;
					ATA[pieceIndex2 * ATAPieceWidth + index][row] = -value;

					xpow *= x;
				}
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
