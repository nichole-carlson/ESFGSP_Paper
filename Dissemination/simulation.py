import numpy as np

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
