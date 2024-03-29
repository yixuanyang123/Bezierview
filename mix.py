import numpy as np
import matplotlib.pyplot as plt

qE = np.zeros((3, 6))  # Initialize qE as a 3x6 zero matrix: 3 rows for 3D points, 6 columns for the weights
# Clearing the current figure
plt.clf()

# -------- functionality
shw = 1
do_prt = 0  # Python plotting: only makes sense when shw > 0
bvout = 0
bbnet = 1

# -------- sizes
cfs = 6  # coefficients per patch (total degree quadratic)
idx1 = np.arange(cfs) # orientation
idx0 = np.array([3, 2, 1, 5, 4, 6])-1
pats = 6  # patches per triangle
n = 32  # evaluation density
dim = 3  # dimension of range (3-space)
tet = 0  # use tetrahedra data

# --------- evaluated BB basis functions
u0 = np.linspace(0, 1, n)
u, v = np.meshgrid(u0, u0)
bbb = [(1 - u - v) ** 2, # not sure if this is correct, but it's not being used
      2 * (1 - u - v) * u,
      u * u,
      2 * (1 - u - v) * v,
      2 * v * u,
      v * v]

# Forms upper-left triangle and makes everything else nan
mask = np.ones((n,n))
for col in range(1,n):
    for row in range(n-col, n):
        mask[row][col] = np.nan

# (3 - 3*c0)/(c0 + 1),   3/(2*(c0 + 1))
#  yields for c0 = -1/2     9, 3
# 9
# 3  1
# 3  1  1
nn = np.array([3, 4])  # allowed valences
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
val = np.array([4, 4, 4, 3, 3]) - 1

# Neighbors
nbr = np.array([
    np.array([1, 4, 3, 5, 4, 1]) - 1,
    np.array([2, 2, 2, 2, 3, 3]) - 1,
    np.array([4, 3, 5, 1, 1, 5]) - 1
]).T

if tet == 1:
    # -------- tet
    V = 3 * np.array([
        [-1, 1, 1, -1],
        [-1, 1, -1, 1],
        [-1, -1, 1, 1]
    ])
    val = np.array([3, 3, 3, 3]) - 1
    nbr = np.array([
        [2, 3, 4],
        [1, 4, 3],
        [4, 1, 2],
        [3, 2, 1]
    ])

nbl = nbr[:, [2, 1, 0]]  # Reorder columns of nbr
dim, vts = V.shape  # Get dimensions of V
fcs, vfc = nbr.shape  # Get dimensions of nbr

bez = np.array([
        [
            [
                [0.0 for i in range(4)] for k in range(cfs)
            ] for l in range(8)
        ] for j in range(pats)
    ])  # Initialize bez as an empty dictionary

fc = np.array([None] * fcs)
ctr = np.zeros((V.shape[0], fcs))
dual = np.array([None] * fcs)

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
    #plt.show()

# --- complete Euclidean part of a single quadratic
if bvout == 1:
    fp = open('mixW.bv', 'w')
for orient in range(2):
    if orient == 0:
        nbs = nbr
        idx = idx1
        if bvout == 1:
            print("group {} odd".format(orient), file=fp)
    else:
        nbs = nbl
        idx = idx0
        if bvout == 1:
            print("group {} even".format(orient), file=fp)

    for ff in range(fcs):
        fc[ff] = V[:, nbs[ff, :]]
        A = fc[ff] @ np.ones((vfc,1))/vfc;
        ctr[:, ff] = A[:,0]
        dual[ff] = np.zeros((V.shape[0], vfc))
        for kk in range(vfc):
            dual[ff][:, kk] = (3 * ctr[:, ff] + fc[ff][:, kk]) / 4

    # --- project vtx neighbors MISSING(not needed for specific geometry)
    for ii in range(vts):
        # project duals -- currently left out since valence 3 or 4
        pass

    for ff in range(fcs):
        for kk in range(vfc):  # k = index inside face if rem(pat,2)==1
            pat = 2 * (kk+1) + (orient+1) - 1
            km = kk - 1
            if km < 0:
                km = kk + vfc - 1
            kp = kk + 1
            if kp >= vfc:
                kp = kk - vfc + 1
            top = nbs[ff, kk]  # top global vertex index
            nxt = nbs[ff, kp]  # bottom left vertex index
            prv = nbs[ff, km]  # bottom right vertex index

            ntp = 0
            for tt in range(fcs):  # neighbor triangle nxt
                for ii in range(vfc):
                    ip = ii + 1
                    if ip >= vfc:
                        ip = 0
                    if nbs[tt, ii] == nxt and nbs[tt, ip] == top:
                        ntp = tt
                        break
                if ntp != 0:
                    break

            ntm = 0
            for jj in range(fcs):  # neighbor triangle prv-top edge
                for nn in range(vfc):
                    np1 = nn + 1
                    if np1 >= vfc:
                        np1 = 0
                    if nbs[jj, nn] == top and nbs[jj, np1] == prv:
                        ntm = jj
                        break
                if ntm != 0:
                    break

            # v10 v00 ff
            # v11 v01
            vf0 = dual[ff][:, kk]
            vfp = dual[ff][:, kp]
            vfm = dual[ff][:, km]
            v0p = dual[ntp][:, ip]
            v1p = dual[ntp][:, ii]
            v0m = dual[ntm][:, nn]
            if pat == 0 and ff == 0:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(vf0[0], vf0[1], vf0[2], color='r', marker='o')
                ax.scatter(v0p[0], v0p[1], v0p[2], color='g', marker='*')
                ax.scatter(vfp[0], vfp[1], vfp[2], color='b', marker='*')
                ax.scatter(v1p[0], v1p[1], v1p[2], color='k', marker='+')
                #plt.show()

            wti_top = wti[val[top - 1] - 2]
            wti_bot = wti[val[nxt - 1] - 2]
            wni_top = wni[val[top - 1] - 2]
            # w_top = 2 * (1 + c0(val(top) - 2)) / 3
            # w_bot = 2 * (1 + c0(val(nxt) - 2)) / 3
            # w_prv = 2 * (1 + c0(val(prv) - 2)) / 3
            w_top = 1
            w_bot = 1
            w_prv = 1
            w1 = np.sqrt(wti_top * wti_bot)
            w5 = w_top
            w4 = w_top * wti_top
            wgt = [w1, w5, (w_top + w_bot + w_prv) / 3,
                   w4, w5,
                   w_top * wni_top]

            # --- assemble by averaging in coeff_i weight_i
            qE[:, 4] = wgt[4] * vf0
            qE[:, 1] = wgt[1] * (vf0 + vfp) / 2
            midopp = wgt[1] * (v0p + v1p) / 2
            qE[:, 3] = wgt[3] * (vf0 + v0p) / 2
            qE[:, 0] = wgt[0] * (vf0 + vfp + v0p + v1p) / 4
            qE[:, 2] = (w_top * vf0 + w_bot * vfp + w_prv * vfm) / 3

            # top is average:  v0p-o + v0m-o = (vf0-o)*2c0
            #  v0p+v0m-vf0(2c0) = 2o-2c0o = 2(1-c0)o
            cc = c0[val[top - 1] - 2]
            tv = (v0p + v0m - 2 * cc * vf0) / (2 * (1 - cc))  # top vertex Euclidian
            qE[:, 5] = wgt[5] * tv
            if pat == 0 and ff == 0 and kk == 0:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(v0m[0], v0m[1], v0m[2], color='c', marker='o')
                ax.scatter(tv[0], tv[1], tv[2], color='r', marker='+')
                #plt.show()

            # --- assemble
            bbase = np.vstack([qE, wgt]).T
            # print(bbase.shape)
            # Show based on condition
            if shw == 4:
                if ff == 0:
                    plt.figure()
                    plt.plot(bbase[:, 0], bbase[:, 1], 'b')
                    #plt.show()

            # print(idx,ff,pat)
            # Export

            bez[ff][pat] = bbase[idx,:][:]
            if bvout == 1:
                fp.write(f"{11}\n{2}\n")
                for cf in range(cfs):  # BB-coeffs of patch
                    fp.write(
                        f"{bez[pat][ff][cf, 0]} {bez[pat][ff][cf, 1]} {bez[pat][ff][cf, 2]} {bez[pat][ff][cf, 3]}\n")
            if shw == 3:
                plt.show(bez[pat][ff], dim, cfs, bbb, mask, 'r')
            if bbnet == 1:
                of = 0.01
                for ii in range(cfs):
                    ww = bbase[ii, 3]
                    xx = bbase[ii, 0] / ww
                    yy = bbase[ii, 1] / ww
                    zz = bbase[ii, 2] / ww
                    ax.text(xx + of, yy + of, zz + of, round(ww,2))
                ids = np.array([[0, 1, 3], [1, 2, 4], [3, 4, 5]])
                for ii in range(3):  # Three subtriangles of bb-net of one quadratic
                    for jj in range(3):  # Each corner of a subtriangle
                        jp = jj + 1
                        if jp == 3:
                            jp = 0

                        ll = [ids[ii][jj], ids[ii][jp]]

                        color = 'k' if ff % 2 == 1 else 'r'
                        ax.plot(bbase[ll, 0] / bbase[ll, 3], bbase[ll, 1] / bbase[ll, 3], bbase[ll, 2] / bbase[ll, 3],
                                 color + '-', linewidth=3)
                        #plt.show()

            if do_prt == 1:
                plt.savefig('tst.png')
plt.show()
