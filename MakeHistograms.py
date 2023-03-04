import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

doHistograms = True
doDFTs = True
doCDFs = True

print("Loading Data")
df = pd.read_csv('out.csv')

columns = df.columns.values.tolist()
columnCount = len(columns)

# ================= Histograms =================
if doHistograms:
    print("Histograms")

    diagramCols = min(int(math.sqrt(columnCount)), 4)
    diagramRows = math.ceil(columnCount / diagramCols)

    fig, ax = plt.subplots(diagramRows, diagramCols, figsize=(15, 20), squeeze=False)

    fig.suptitle('Histograms', fontsize=16)

    for i in range(columnCount):
        print("  " + columns[i])
        ax[i%diagramRows, math.floor(i/diagramRows)].hist(df[columns[i]], bins=100)
        ax[i%diagramRows, math.floor(i/diagramRows)].set_title(columns[i])

    plt.tight_layout()
    fig.savefig("_histograms.png", bbox_inches='tight')
    fig.savefig("_histograms.pdf", bbox_inches='tight')

# ================= DFTs =================
if doDFTs:
    print("DFT")

    graphsPerCell = 4

    columns = df.columns.values.tolist()
    columnCount = len(columns)
    graphCount = int(columnCount / graphsPerCell)

    diagramCols = min(int(math.sqrt(graphCount)), 4)
    diagramRows = math.ceil(graphCount / diagramCols)    

    fig, ax = plt.subplots(diagramRows, diagramCols, figsize=(15, 10), squeeze=False)

    fig.suptitle('DFTs', fontsize=16)

    numDFTsAvgd = 1000

    for i in range(graphCount):
        print("  " + columns[i*graphsPerCell])

        ax[i%diagramRows, math.floor(i/diagramRows)].set_title(columns[i*graphsPerCell])

        lineStyles = ['-',':',':', ':']

        for j in range(graphsPerCell):
            data = df[columns[i*graphsPerCell + j]].to_numpy()

            AvgDFT = None
            for dftIndex in range(numDFTsAvgd):
                startIndex = int(dftIndex * len(data) / numDFTsAvgd)
                stopIndex = int((dftIndex+1) * len(data) / numDFTsAvgd)
                stopIndex = min(stopIndex, len(data))
                DFT = np.fft.fft(data[startIndex:stopIndex])

                # zero out DC and only show the positive frequencies
                DFT[0] = 0
                DFT = DFT[1:int(math.ceil(len(DFT)/2))]

                #DFT = np.fft.fftshift(DFT)
                DFT = np.log(1 + abs(DFT))
                if dftIndex == 0:
                    AvgDFT = DFT
                else:
                    alpha = 1.0 / float(dftIndex + 1)
                    AvgDFT = AvgDFT * (1.0 - alpha) + DFT * alpha

            line, = ax[i%diagramRows, math.floor(i/diagramRows)].plot(AvgDFT, lineStyles[j])
            line.set_label(columns[i*graphsPerCell+j])

        labelsPos = []
        labels = []
        labelCount = 11
        xwidth = len(data) / numDFTsAvgd
        for labelIndex in range(labelCount):
            labelsPos.append(0.5 * xwidth * labelIndex / (labelCount-1))
            labels.append(labelIndex / (labelCount-1))

        ax[i%diagramRows, math.floor(i/diagramRows)].set_xticks(labelsPos, labels=labels)
        ax[i%diagramRows, math.floor(i/diagramRows)].legend()

    plt.tight_layout()
    fig.savefig("_DFTs.png", bbox_inches='tight')
    fig.savefig("_DFTs.pdf", bbox_inches='tight')

# ================= CDFs =================
if doCDFs:
    print("CDFs")

    df = pd.read_csv('cdf.csv')

    graphsPerCell = 3

    columns = df.columns.values.tolist()
    columnCount = len(columns)
    graphCount = int(columnCount / graphsPerCell)

    diagramCols = min(int(math.sqrt(graphCount)), 4)
    diagramRows = math.ceil(graphCount / diagramCols)

    fig, ax = plt.subplots(diagramRows, diagramCols, figsize=(15, 10), squeeze=False)

    fig.suptitle('CDFs', fontsize=16)

    for i in range(graphCount):
        print("  " + columns[i*graphsPerCell])
        ax[i%diagramRows, math.floor(i/diagramRows)].set_title(columns[i*graphsPerCell])
        for j in reversed(range(graphsPerCell)):
            line, = ax[i%diagramRows, math.floor(i/diagramRows)].plot(df[columns[i*graphsPerCell + j]])
            line.set_label(columns[i*graphsPerCell+j])
        ax[i%diagramRows, math.floor(i/diagramRows)].legend()

    plt.tight_layout()
    fig.savefig("_cdf.png", bbox_inches='tight')
    fig.savefig("_cdf.pdf", bbox_inches='tight')

    # Do the CDF error graphs
    print("CDF Errors")

    fig, ax = plt.subplots(diagramRows, diagramCols, figsize=(15, 10), squeeze=False)

    fig.suptitle('CDF Errors', fontsize=16)

    for i in range(graphCount):
        print("  " + columns[i*graphsPerCell])
        ax[i%diagramRows, math.floor(i/diagramRows)].set_title(columns[i*graphsPerCell] + " Error")
        for j in reversed(range(1, graphsPerCell)):
            line, = ax[i%diagramRows, math.floor(i/diagramRows)].plot(df[columns[i*graphsPerCell + j]] - df[columns[i*graphsPerCell + 0]])
            line.set_label(columns[i*graphsPerCell+j])
        ax[i%diagramRows, math.floor(i/diagramRows)].legend()

    plt.tight_layout()
    fig.savefig("_cdferror.png", bbox_inches='tight')
    fig.savefig("_cdferror.pdf", bbox_inches='tight')
