import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

n_a = 1000
n_b = 1000
n_image = n_a + n_b 
n_pixel = 256
n_iter = 100
image_size = int(np.sqrt(n_pixel))
center_size = 8


# Simulate data --------------------------------------------------
# 1. Generate group indicator, 1000 in group A; 1000 in group B 
# 2. Generate the beta as 16*16 matrix, it has 1 in the center 8*8 and 0 in 
#    other areas. 
# 3. Generate the exponential correlation matrix 
# 4. Use the corr matrix to generate the epsilon as a multivariate normal. 
#    Epsilon should be (100, 2000, 256)
# 5. Broadcast and generate y 

group_ind = np.concatenate((np.ones(n_a), np.zeros(n_b)))
group_ind = group_ind[np.newaxis, :, np.newaxis]

beta = np.zeros((image_size, image_size))
central_start = (image_size - center_size) // 2
central_end = central_start + center_size
beta[central_start:central_end, central_start:central_end] = 1
beta = beta.flatten()

def exp_corr_mat(n, rate):
    indices = np.arange(n)
    dist_mat = np.abs(indices - indices[:, np.newaxis])
    corr_mat = np.exp(-rate * dist_mat)
    return corr_mat

corr_mat = exp_corr_mat(n_pixel, rate=1)

epsilon = np.random.multivariate_normal(
    mean=np.zeros(n_pixel), cov=corr_mat, size=(n_iter, n_image)
)
y = beta * group_ind + epsilon


# Visualization --------------------------------------------------------
plt.imshow(y[0, 768, :].reshape((16, 16)), cmap="viridis")
plt.colorbar(label="Value")
plt.show()



# VBM ------------------------------------------------------------------
p_values = np.zeros((n_iter, n_pixel))

for iter_index in range(n_iter):
    for pixel_index in range(n_pixel):
        y_pixel = y[iter_index, :, pixel_index]
        model = sm.GLM(y_pixel, group_ind[0, :, :], family = sm.families.Gaussian())
        p_values[iter_index, pixel_index] = model.fit().pvalues[0]
    p_values[iter_index, :] = multipletests(p_values[iter_index, :], method="bonferroni")[1]

# calculate the perc of p-values < 0.05 for each pixel
pvals_summ = np.sum(p_values < 0.05, axis=0) / n_iter * 100

plt.imshow(
    pvals_summ.reshape((image_size, image_size)), 
    cmap="gray_r", vmin=0, vmax =100, interpolation="none"
)
plt.title("Corrected p-values")
plt.colorbar(label='Percentage of p-values < 0.05')
plt.xlabel("Pixel X-coordinate")
plt.ylabel("Pixel Y-coordinate")
plt.show()
