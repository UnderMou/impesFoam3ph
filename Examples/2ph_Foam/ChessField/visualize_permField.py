import cv2
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("Agg")

min_val, max_val = 1e-13, 4e-13

# Step 1: Load the image
image_path = "perm_quad.png"  # Replace with the path to your image
image = cv2.imread(image_path)

# Step 2: Convert the image to grayscale
gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# Step 3: The grayscale image is already a NumPy array, but let’s make it explicit
gray_matrix = np.array(gray_image)
print(gray_matrix.shape)
normalized_gray = gray_matrix / 255.0  # Values from 0 (black) to 1 (white)
real_value_matrix = normalized_gray * (max_val - min_val) + min_val

plt.figure(figsize=(6, 6))
plt.imshow(real_value_matrix, cmap='gray')  # Display the matrix as an image in grayscale
plt.axis("off")  # Hide axes for a cleaner look
plt.savefig('scanned_perm_quad.png', dpi=300)

plt.figure(figsize=(6, 6))
im = plt.imshow(real_value_matrix, cmap="gray", vmin=min_val, vmax=max_val)  # Use real value range for color bar
plt.axis("off")
cbar = plt.colorbar(im, fraction=0.04, pad=0.04, shrink=0.7, aspect=10)
cbar.set_label(r"K [$m^2$]")
cbar.ax.tick_params(labelsize=8)
plt.savefig('scanned_perm_quad.png', dpi=300)

# K0 = real_value_matrix[0,:]
# print(K0)
# x0 = np.linspace(0,1,len(K0))
# plt.plot(x0,K0)
# plt.savefig('K0.png',dpi=300)
