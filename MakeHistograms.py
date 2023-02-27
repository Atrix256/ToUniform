import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

df = pd.read_csv('out.csv')

columns = df.columns.values.tolist()
columnCount = len(columns)

rows = int(columnCount / 4)

# ================= Histograms =================

print("Histograms")

fig, ax = plt.subplots(rows, math.ceil(columnCount/rows), figsize=(15, 10))

for i in range(columnCount):
    print("  " + columns[i])
    ax[i%rows, math.floor(i/rows)].hist(df[columns[i]], bins=100)
    ax[i%rows, math.floor(i/rows)].set_title(columns[i])

plt.tight_layout()
fig.savefig("_histograms.png", bbox_inches='tight')

# ================= DFTs =================

print("DFT")

fig, ax = plt.subplots(rows, math.ceil(columnCount/rows), figsize=(15, 10))

numDFTsAvgd = 1000

for i in range(columnCount):
    print("  " + columns[i])
    data = df[columns[i]].to_numpy()

    AvgDFT = None
    for dftIndex in range(numDFTsAvgd):
        startIndex = int(dftIndex * len(data) / numDFTsAvgd)
        stopIndex = int((dftIndex+1) * len(data) / numDFTsAvgd)
        stopIndex = min(stopIndex, len(data))
        DFT = np.fft.fft(data[startIndex:stopIndex])
        DFT[0] = 0 # remove DC
        DFT = np.fft.fftshift(DFT)
        DFT = np.log(1 + abs(DFT))
        if dftIndex == 0:
            AvgDFT = DFT
        else:
            alpha = 1.0 / float(dftIndex + 1)
            AvgDFT = AvgDFT * (1.0 - alpha) + DFT * alpha

    ax[i%rows, math.floor(i/rows)].plot(AvgDFT)
    ax[i%rows, math.floor(i/rows)].set_title(columns[i])
    #ax[i%rows, math.floor(i/rows)].get_xaxis().set_visible(False)
    #ax[i%rows, math.floor(i/rows)].get_yaxis().set_visible(False)
    #freq = np.fft.fftfreq(DFT.size, d=0.1)
    #freq = np.fft.fftshift(freq)
    #print(data.size)
    #print(freq.size)
    #ax[i%rows, math.floor(i/rows)].set_xticks(freq)

plt.tight_layout()
fig.savefig("_DFTs.png", bbox_inches='tight')

# ================= CDFs =================

print("CDFs")

df = pd.read_csv('cdf.csv')

graphsPerCell = 2

columns = df.columns.values.tolist()
columnCount = len(columns)
graphCount = int(columnCount / graphsPerCell)

rows = int(graphCount / 4)

fig, ax = plt.subplots(rows, math.ceil(graphCount/rows), figsize=(15, 10))

for i in range(graphCount):
    print("  " + columns[i*graphsPerCell])
    ax[i%rows, math.floor(i/rows)].set_title(columns[i*graphsPerCell])
    for j in reversed(range(graphsPerCell)):
        line, = ax[i%rows, math.floor(i/rows)].plot(df[columns[i*graphsPerCell + j]])
        line.set_label(columns[i*graphsPerCell+j])
    ax[i%rows, math.floor(i/rows)].legend()

plt.tight_layout()
fig.savefig("_cdf.png", bbox_inches='tight')
