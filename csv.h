#pragma once

#include <string>
#include <vector>

struct Column
{
	std::string label;
	std::vector<float> values;
};
typedef std::vector<Column> CSV;

inline void WriteCSV(const CSV& csv, const char* fileName)
{
	FILE* file = nullptr;
	fopen_s(&file, fileName, "wb");

	if (csv.size() > 0)
	{
		// write header
		size_t maxColSize = 0;
		bool first = true;
		for (const Column& column : csv)
		{
			fprintf(file, "%s\"%s\"", first ? "" : ",", column.label.c_str());
			maxColSize = std::max(maxColSize, column.values.size());
			first = false;
		}
		fprintf(file, "\n");

		// write data
		for (size_t i = 0; i < maxColSize; ++i)
		{
			bool first = true;
			for (const Column& column : csv)
			{
				if (column.values.size() > i)
					fprintf(file, "%s\"%f\"", first ? "" : ",", column.values[i]);
				else
					fprintf(file, "%s\"\"", first ? "" : ",");
				first = false;
			}
			fprintf(file, "\n");
		}
	}

	fclose(file);
}