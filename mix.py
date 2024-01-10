import numpy as np
import matplotlib.pyplot as plt

# Clearing the current figure
plt.clf()

# -------- functionality
shw = 1
do_prt = 0  # Python plotting: only makes sense when shw > 0
bvout = 1
bbnet = 1

# -------- sizes
cfs = 6  # coefficients per patch (total degree quadratic)
idx1 = np.arange(1, cfs + 1)  # orientation
idx0 = [3, 2, 1, 5, 4, 6]
pats = 6  # patches per triangle
n = 32  # evaluation density
dim = 3  # dimension of range (3-space)
tet = 0  # use tetrahedra data

# --------- evaluated BB basis functions
u0 = np.linspace(0, 1, n)
u, v = np.meshgrid(u0, u0)
bbb = {
    1: (1 - u - v)**2,
    2: 2 * (1 - u - v) * u,
    3: u * u,
    4: 2 * (1 - u - v) * v,
    5: 2 * v * u,
    6: v * v
}
mask = np.ones((n, n))  # suppress half of the 4-sided patch
bidx = np.triu_indices(n, k=1)
mask[bidx] = np.nan

# Valences
nn = np.array([3, 4])
c0 = np.cos(2 * np.pi / nn)

# Rational weights
wni = 3 * (1 - c0) / (c0 + 1)
wti = 3 / (2 * (c0 + 1))

# -------- GEOMETRY+CONNECTIVITY:  double simplex
# Vertices
V = np.array([
    np.cos(2 * np.pi * np.array([0, 1, 2]) / 3).tolist() + [0, 0],
    np.sin(2 * np.pi * np.array([0, 1, 2]) / 3).tolist() + [0, 0],
    [0, 0, 0, -1, 1]
])
val = np.array([4, 4, 4, 3, 3])

# Neighbors
nbr = np.array([
    [1, 4, 3, 5, 4, 1],
    [2, 2, 2, 2, 3, 3],
    [4, 3, 5, 1, 1, 5]
]).T

if tet == 1:
    # -------- tet
    V = 3 * np.array([
        [-1, 1, 1, -1],
        [-1, 1, -1, 1],
        [-1, -1, 1, 1]
    ])
    val = np.array([3, 3, 3, 3])
    nbr = np.array([
        [2, 3, 4],
        [1, 4, 3],
        [4, 1, 2],
        [3, 2, 1]
    ])

nbl = nbr[:, [2, 1, 0]]  # Reorder columns of nbr
dim, vts = V.shape  # Get dimensions of V
fcs, vfc = nbr.shape  # Get dimensions of nbr

# ---- draw funnel
if shw > 0:
    clr = ['y', 'c', 'r']  # Colors
    fidx = [1, 2, 5]  # Indices for funnel drawing

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for jj in range(3):
        ii = fidx[jj] - 1  # 0-based indexing

    ax.view_init(elev=-V[0, -1], azim=V[1, -1])  # 0-based indexing
    ax.axis('equal')
    plt.show()

# --- complete Euclidean part of a single quadratic
if bvout == 1:
    fp = open('mixW.bv', 'w')
