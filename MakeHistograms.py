import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

df = pd.read_csv('out.csv')
rows = 5

columns = df.columns.values.tolist()
columnCount = len(columns)

# ================= Histograms =================

print("Histograms")

fig, ax = plt.subplots(rows, math.ceil(columnCount/rows), figsize=(15, 10))

for i in range(columnCount):
    print(columns[i])
    ax[i%rows, math.floor(i/rows)].hist(df[columns[i]], bins=100)
    ax[i%rows, math.floor(i/rows)].set_title(columns[i])

plt.tight_layout()
fig.savefig("histograms.png", bbox_inches='tight')

# ================= DFTs =================

print("DFT")

fig, ax = plt.subplots(rows, math.ceil(columnCount/rows), figsize=(15, 10))

for i in range(columnCount):
    print(columns[i])
    data = df[columns[i]].to_numpy()
    DFT = np.fft.fft(data)
    DFT[0] = 0 # remove DC
    DFT = np.fft.fftshift(DFT)
    DFT = np.log(1 + abs(DFT))
    ax[i%rows, math.floor(i/rows)].plot(DFT)
    ax[i%rows, math.floor(i/rows)].set_title(columns[i])
    #ax[i%rows, math.floor(i/rows)].get_xaxis().set_visible(False)
    #ax[i%rows, math.floor(i/rows)].get_yaxis().set_visible(False)
    #freq = np.fft.fftfreq(DFT.size, d=0.1)
    #freq = np.fft.fftshift(freq)
    #print(data.size)
    #print(freq.size)
    #ax[i%rows, math.floor(i/rows)].set_xticks(freq)

plt.tight_layout()
fig.savefig("DFTs.png", bbox_inches='tight')