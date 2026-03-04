import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# 1. Generate some sample discrete data (integers from 0 to 4)
# data = np.random.randint(0, 5, size=(10, 10))
data = np.random.randint(0, 4, size=(5, 5))
print(data)

# 2. Define the discrete colormap and the normalization
# List the specific colors you want to use
# colors = ['navy', 'crimson', 'limegreen', 'gold', 'purple']
colors = ['navy', 'crimson', 'limegreen', 'gold']
cmap = mcolors.ListedColormap(colors)
# Define the boundaries for each color.
# There should be N+1 boundaries for N colors.
# bounds = [0, 1, 2, 3, 4, 5]
bounds = [0, 1, 2, 3, 4]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# 3. Create the plot and display the data
fig, ax = plt.subplots(figsize=(7, 4))
# Use imshow for a 2D array plot
image = ax.imshow(data, cmap=cmap, norm=norm)

# 4. Add the horizontal colorbar
# Specify 'horizontal' orientation and define explicit ticks
cbar = fig.colorbar(image, ax=ax, orientation='horizontal', ticks=bounds[:-1], pad=0.1)
cbar.set_label('Discrete Values')

# Set the tick labels to be centered within the color segments
print("bounds=",bounds)
cbar.set_ticks([bound + 0.5 for bound in bounds[:-1]])
# cbar.set_ticklabels(['0', '1', '2', '3', '4'])
cbar.set_ticklabels(['0', '1', '2', '3'])

# 5. Display the plot
plt.show()