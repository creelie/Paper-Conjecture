"""Figures for 'Degrees of Prolongation Varieties of Monomial Planar Systems'."""
import sys
from fractions import Fraction
from math import gcd
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

OUT = sys.argv[1] if len(sys.argv) > 1 else "figures"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})

def jets(p, q, L):
    """Supports of A_l, B_l as dicts {(i,j): coeff} (i = a-exponent, j = b-exponent)."""
    def D(h):
        out = {}
        for (i, j), c in h.items():
            if i: out[(i-1, j+p)] = out.get((i-1, j+p), 0) + c*i
            if j: out[(i+q, j-1)] = out.get((i+q, j-1), 0) + c*j
        return out
    A, B = [{(1, 0): 1}], [{(0, 1): 1}]
    for _ in range(L):
        A.append(D(A[-1])); B.append(D(B[-1]))
    return A, B

# Figure 1: defect bars for (p,q)=(3,4)
def fig_defect():
    p, q = 3, 4; c = p*q-1
    fig = plt.figure(figsize=(7.2, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("viridis")
    vals = {}
    for r in range(p+1):
        for s in range(q+1):
            vals[(r, s)] = c*(Fraction(r*(p-r), p+1) + Fraction(s*(q-s), q+1))
    vmax = float(max(vals.values()))
    for (r, s), v in sorted(vals.items(), key=lambda t: -(t[0][0]+t[0][1])):
        v = float(v)
        col = cmap(0.15 + 0.75*v/vmax) if v > 0 else (0.8, 0.8, 0.8, 1)
        ax.bar3d(r-0.32, s-0.32, 0, 0.64, 0.64, max(v, 0.02), color=col, edgecolor="k", linewidth=0.4, shade=True)
    ax.set_xlabel(r"$r \equiv m \ (\mathrm{mod}\ 4)$", labelpad=6)
    ax.set_ylabel(r"$s \equiv m \ (\mathrm{mod}\ 5)$", labelpad=6)
    ax.set_zlabel("loss", labelpad=4)
    ax.set_xticks(range(p+1)); ax.set_yticks(range(q+1))
    ax.view_init(elev=26, azim=-58)
    ax.set_box_aspect((1, 1.25, 0.8))
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_defect.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

# Figure 2: interlacing ribbons for (p,q)=(2,3)
def fig_interlace():
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    p, q = 2, 3; N, M = p+1, q+1
    u, v = sp.symbols("u v")
    X = lambda h: sp.expand(sp.Rational(1, N)*v**q*sp.diff(h, u) - sp.Rational(1, M)*u**p*sp.diff(h, v))
    phi, psi = v, u
    xs = np.linspace(0, 1, 801)
    U = xs**(1/N); V = (1-xs)**(1/M)
    fig = plt.figure(figsize=(8.4, 5.6))
    ax = fig.add_subplot(111, projection="3d")
    orders = range(3, 10)
    for _ in range(min(orders)):
        phi, psi = X(phi), X(psi)
    for l in orders:
        rA, sA = (-l) % N, (1-l) % M
        rB, sB = (1-l) % N, (-l) % M
        for h, col, off, (rr, ss) in ((phi, "#1f5fbf", 0.0, (rA, sA)), (psi, "#c0392b", 0.0, (rB, sB))):
            # F_l = phi_l / (u^rho v^sigma) is a polynomial in x with the same zeros in (0,1)
            hh = sp.cancel(h/(u**rr*v**ss))
            f = sp.lambdify((u, v), hh, "numpy")
            y = np.asarray(f(U, V), dtype=float) * np.ones_like(xs)
            y = y/np.max(np.abs(y))
            yy = np.full_like(xs, l+off)
            verts = [list(zip(xs, yy, y)) + [(1, l+off, 0), (0, l+off, 0)]]
            ax.plot(xs, yy, y, color=col, lw=1.4)
            if col == "#1f5fbf":
                ax.add_collection3d(Poly3DCollection([[(0, l-0.35, 0), (1, l-0.35, 0), (1, l+0.35, 0), (0, l+0.35, 0)]], facecolor="0.85", alpha=0.35, edgecolor="none"))
                ax.plot([0, 1], [l, l], [0, 0], color="0.4", lw=0.6)
            sgn = np.sign(y)
            idx = np.where((sgn[:-1]*sgn[1:] < 0))[0]
            zx = [xs[i] - y[i]*(xs[i+1]-xs[i])/(y[i+1]-y[i]) for i in idx]
            zx = [z for z in zx if 1e-3 < z < 1-1e-3]
            ax.scatter(zx, [l+off]*len(zx), [0]*len(zx), color=col, s=20, depthshade=False, edgecolor="k", linewidth=0.4)
        phi, psi = X(phi), X(psi)
    ax.set_xlabel(r"$x=u^{3}$", labelpad=8); ax.set_ylabel(r"order $l$", labelpad=8)
    ax.set_zlabel("normalised value", labelpad=6)
    ax.set_xlim(0, 1); ax.set_ylim(2.6, 9.4); ax.set_zlim(-1.6, 1.6)
    ax.set_yticks(list(orders)); ax.set_yticklabels([str(l) for l in orders])
    ax.set_zticks([-1, 0, 1]); ax.set_xticks([0, 0.5, 1])
    ax.view_init(elev=38, azim=-72)
    ax.set_box_aspect((1.5, 2.2, 0.55))
    ax.plot([], [], color="#1f5fbf", lw=1.4, label=r"$F_l$ (zeros of $A_l$)")
    ax.plot([], [], color="#c0392b", lw=1.4, label=r"$G_l$ (zeros of $B_l$)")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_interlace.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

# Figure 3: Newton diagram at [0:1:0] for (p,q)=(2,3), m=4
def fig_newton():
    p, q, m = 2, 3, 4; l = m+1
    g = gcd(p+1, q+1); al, be, de = (p+1)//g, (q+1)//g, (p*q-1)//g
    A, B = jets(p, q, l)
    WA, WB = al+l*de, be+l*de
    fpts, gpts = set(), set()
    for W, top, pts in ((WA, A[l], fpts), (WB, B[l], gpts)):
        for (i, j) in top: pts.add((0, j))
        for t in range(0, l):
            for poly, w0 in ((A[t], al+t*de), (B[t], be+t*de)):
                for (i, j) in poly: pts.add((W-w0, j))
        pts.add((W, 0))
    lam = l % (p+1)
    fP, eP, fQ, eQ = p+1-lam, al-be+(lam-1)*de, p+2-lam, (lam-1)*de
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    E, J = 13, 6
    for (e, j) in fpts:
        if e <= E and j <= J: ax.plot(e, j, "o", mfc="none", mec="#1f5fbf", ms=9, mew=1.4)
    for (e, j) in gpts:
        if e <= E and j <= J: ax.plot(e, j, "s", mfc="#c0392b", mec="#c0392b", ms=4.5, alpha=0.85)
    ax.fill([0, eP, E+0.5, E+0.5, 0], [fP, 0, 0, J+0.5, J+0.5], color="#1f5fbf", alpha=0.07, lw=0)
    ax.plot([0, eP], [fP, 0], color="#1f5fbf", lw=2, label=r"edge of $f$: $(0,%d)$ to $(%d,0)$" % (fP, eP))
    ax.plot([0, eQ], [fQ, 0], color="#c0392b", lw=2, ls="--", label=r"edge of $g$: $(0,%d)$ to $(%d,0)$" % (fQ, eQ))
    ax.set_xlim(-0.5, E+0.5); ax.set_ylim(-0.5, J+0.5)
    ax.set_xlabel(r"exponent $e$ of $Z$"); ax.set_ylabel(r"exponent $j$ of $b$")
    ax.set_xticks(range(0, E+1)); ax.set_yticks(range(0, J+1))
    ax.grid(color="0.88", lw=0.6); ax.set_axisbelow(True)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_newton.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    return fP, eP, fQ, eQ

fig_defect(); fig_interlace(); print(fig_newton())
