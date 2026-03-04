import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

# Create a 5x5 grid of subplots
fig, axes = plt.subplots(5, 5, figsize=(15, 15), constrained_layout=True)

# Flatten the 2D array of axes to 1D for easy iteration
axes_flat = axes.flatten()

# Dummy loop: Replace this with your actual image loading logic
# Assuming you have a list of image paths: image_files = ["img1.png", "img2.png", ...]
for i, ax in enumerate(axes_flat):
    # --- Example: Create dummy images (replace with: img = mpimg.imread(image_files[i])) ---
    dummy_img = np.random.rand(100, 100)
    ax.imshow(dummy_img, cmap='gray')
    # ----------------------------------------------------------------------------------------

    # Add axis labels and title
    ax.set_title(f'Image {i+1}')
    ax.set_xlabel('X-pixel')
    ax.set_ylabel('Y-pixel')
    
    # Optional: Turn off ticks if you only want the label
    # ax.set_xticks([])
    # ax.set_yticks([])

plt.show()
