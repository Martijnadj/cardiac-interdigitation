import matplotlib.pyplot as plt
import numpy as np
import csv
import seaborn as sns





mask_file = '../../../parameters/double_masks/files/mask29.dat'

with open(mask_file, newline='') as csvfile:
    data_mask = list(csv.reader(csvfile))

x_max = 460
y_max = 310

mask_array = np.zeros((x_max, y_max))
for i in range (0,len(data_mask)):
    x = data_mask[i][0]
    y = data_mask[i][1]
    mask_array[int(x)][int(y)] = True



corr = np.random.rand(460, 310)

print(np.shape(corr))
print(np.shape(mask_array))

mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    ax = sns.heatmap(corr.T, mask=(1-mask_array.T), vmax=.3, square=True,  cmap="YlGnBu")
    plt.show()







with open('SF_file.txt', newline='') as csvfile:
    data_SF = list(csv.reader(csvfile))

mask_file = '../../../parameters/double_masks/files/mask29.dat'
with open(mask_file, newline='') as csvfile:
    data_mask = list(csv.reader(csvfile))


x_max = 460
y_max = 310

mask_array = np.zeros((x_max, y_max))
for i in range (0,len(data_mask)):
    x = data_mask[i][0]
    y = data_mask[i][1]
    value = data_mask[i][2]
    mask_array[int(x)][int(y)] = bool(value)

SF_array = np.zeros((x_max, y_max))
for i in range (0,len(data_SF)):
    x = data_SF[i][0]
    y = data_SF[i][1]
    value = data_SF[i][2]
    SF_array[int(x)][int(y)] = value


SF_plot = np.zeros((x_max, y_max))



for x in range (0,x_max):
    for y in range (0,y_max):
        if(mask_array[x][y]):
            SF_plot[x][y] = SF_array[x][y]
        else:
            SF_plot[x][y] = np.nan

ax = sns.heatmap(SF_plot.T)
plt.show()

#print(data_SF)




