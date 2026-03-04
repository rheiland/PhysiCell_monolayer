import matplotlib.pyplot as plt
import numpy as np

# 1. Create dummy image data (replace with your PNG loading logic)
# In reality, you'd load images using: image = plt.imread('path.png')
rows, cols = 5, 5
fig, axs = plt.subplots(rows, cols, figsize=(10, 10), sharex=True, sharey=True)

# 2. Loop through axes and display images
for i in range(rows):
    for j in range(cols):
        # Generate dummy data for demonstration
        dummy_img = np.random.rand(100, 100) 
        
        # Plot the image
        axs[i, j].imshow(dummy_img, cmap='gray')
        
        # Remove tick labels for inner plots, keep on outer edges
        if i == rows - 1:
            axs[i, j].set_xlabel(f'Col {j+1}')
        if j == 0:
            axs[i, j].set_ylabel(f'Row {i+1}')
        
        # Hide tick marks for cleaner layout
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])

# 3. Adjust spacing between subplots
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.suptitle('5x5 Image Grid', fontsize=16)
plt.show()

