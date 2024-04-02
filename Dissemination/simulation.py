import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

na = 1000  # number of group A
nb = 1000  # number of group B
s = 256  # number of pixels

# set up z
z = np.concatenate([np.ones(na), np.zeros(nb)])

# set up beta
beta = np.zeros((16, 16))
central_start = (16 - 8) // 2  # Calculate the starting index
central_end = central_start + 8  # Calculate the ending index
beta[central_start:central_end, central_start:central_end] = 1
beta = beta.flatten()


# set up correlation matrix
def exp_corr_mat(n, rate):
    indices = np.arange(n)
    dist_mat = np.abs(indices - indices[:, np.newaxis])
    corr_mat = np.exp(-rate * dist_mat)
    return corr_mat


corr_mat = exp_corr_mat(s, rate=1)

# set up y
y = np.zeros((na + nb, s))
for i in range(na + nb):
    y[i, :] = np.random.multivariate_normal(z[i] * beta, corr_mat)

# VBM
X = z[:, np.newaxis]  # Makes z a 2D array with a single column

results = []
for pixel in range(y.shape[1]):
    # The response variable for the current pixel across all images
    y_pixel = y[:, pixel]

    model = sm.GLM(y_pixel, X, family=sm.families.Gaussian())
    res = model.fit()
    results.append(res)

vbm_pvals = np.array([res.pvalues[0] for res in results])
vbm_coefs = np.array([res.params[0] for res in results])

# correct p-values for multiple comparison
vbm_correct_pvals = multipletests(vbm_pvals, method="bonferroni")[1]

# coefficients is a 1D numpy array with the estimated parameters
# We reshape it to the 2D format of the images (16x16 in your case)
vbm_coefs_2d = vbm_coefs.reshape((16, 16))

# plot the coefs, black for 1, white for 0
cmap = plt.cm.gray_r
plt.imshow(vbm_coefs_2d, cmap=cmap, interpolation="none")
plt.colorbar(label="Estimated Parameter Value")
plt.title("2D Visualization of Estimated Parameters")
plt.show()

# plot for corrected p-values, black for < 0.05
vbm_correct_pvals_2d = vbm_correct_pvals.reshape((16, 16))
plt.imshow(vbm_correct_pvals_2d, cmap="gray", interpolation="none")
plt.title("Corrected p-values")
plt.xlabel("Pixel X-coordinate")
plt.ylabel("Pixel Y-coordinate")
plt.show()
