import numpy as np
from scipy.optimize import fsolve

def hat(y):
    return np.array([[0, -y[2], y[1]], [y[2], 0, -y[0]], [-y[1], y[0], 0]])

def SingleSectionCR(Force, length, resolution):
    global p, R, j, n, m, v, u, q, w, vs, us, vt, ut, qt, wt, vst, ust, vh, uh, vsh, ush, qh, wh, nLL, mLL, x, y, z, X, Y, Z

    # Parameters
    force = Force
    L = length
    N = resolution
    E = 207e6
    r = 0.01
    rt1 = np.array([0.01, 0, 0])
    rt2 = np.array([0, 0.01, 0])
    rho = 8000
    g = np.array([0, 0, 0])
    Bse = np.zeros((3, 3))
    Bbt = 1e-6 * np.eye(3)
    C = 0.03 * np.eye(3)
    dt = 0.015
    alpha = -0.2
    STEPS = 300
    vstar = lambda s: np.array([0, 0, 1])
    ustar = lambda s: np.array([0, 0, 0])
    vsstar = lambda s: np.array([0, 0, 1])
    usstar = lambda s: np.array([0, 0, 0])

    # Boundary Conditions
    p = {i: np.array([0, 0, 0]) for i in range(STEPS)}
    R = {i: np.eye(3) for i in range(STEPS)}
    q = {i: np.array([0, 0, 0]) for i in range(STEPS)}
    w = {i: np.array([0, 0, 0]) for i in range(STEPS)}

    nL = np.zeros(3)
    mL = np.zeros(3)

    # Dependent Parameter Calculations
    A = np.pi * r**2
    J = np.diag([np.pi * r**4 / 4, np.pi * r**4 / 4, np.pi * r**4 / 2])
    G = E / (2 * (1 + 0.3))
    Kse = np.diag([G * A, G * A, E * A])
    Kbt = np.diag([E * J[0, 0], E * J[1, 1], G * J[2, 2]])
    ds = L / (N - 1)
    c0 = (1.5 + alpha) / (dt * (1 + alpha))
    c1 = -2 / dt
    c2 = (0.5 + alpha) / (dt * (1 + alpha))
    d1 = alpha / (1 + alpha)

    # Main Simulation
    i = 1
    fsolve(staticIVP, np.zeros(6))
    applyStaticBDFalpha()

    for i in range(2, STEPS):
        if i < 5:
            Tt1 = 0
            Tt2 = 0
        else:
            print("Force")
            Tt1 = force
            Tt2 = 0

        fsolve(dynamicIVP, np.concatenate((n[i-1], m[i-1])))
        applyDynamicBDFalpha()

    visualize1()

def applyStaticBDFalpha():
    global vh, uh, vsh, ush, qh, wh, v, u, vs, us, q, w
    for j in range(N-1):
        vh[i+1, j] = (c1 + c2) * v[i, j]
        uh[i+1, j] = (c1 + c2) * u[i, j]
        vsh[i+1, j] = (c1 + c2) * vs[i, j]
        ush[i+1, j] = (c1 + c2) * us[i, j]
        qh[i+1, j] = np.zeros(3)
        wh[i+1, j] = np.zeros(3)
        q[i, j] = np.zeros(3)
        w[i, j] = np.zeros(3)

def applyDynamicBDFalpha():
    global vh, uh, vsh, ush, qh, wh, v, u, vs, us, q, w, vt, ut, qt, wt, vst, ust
    for j in range(N-1):
        vh[i+1, j] = c1 * v[i, j] + c2 * v[i-1, j] + d1 * vt[i, j]
        uh[i+1, j] = c1 * u[i, j] + c2 * u[i-1, j] + d1 * ut[i, j]
        vsh[i+1, j] = c1 * vs[i, j] + c2 * vs[i-1, j] + d1 * vst[i, j]
        ush[i+1, j] = c1 * us[i, j] + c2 * us[i-1, j] + d1 * ust[i, j]
        qh[i+1, j] = c1 * q[i, j] + c2 * q[i-1, j] + d1 * qt[i, j]
        wh[i+1, j] = c1 * w[i, j] + c2 * w[i-1, j] + d1 * wt[i, j]

def staticIVP(G):
    global n, m, p, R, vs, us, v, u
    n[i] = G[:3]
    m[i] = G[3:]

    for j in range(N-1):
        ps, Rs, ns, ms, us[i, j], vs[i, j], v[i, j], u[i, j] = staticODE(p[i, j], R[i, j], n[i, j], m[i, j])
        p[i, j+1] = p[i, j] + ds * ps
        R[i, j+1] = R[i, j] + ds * Rs
        n[i, j+1] = n[i, j] + ds * ns
        m[i, j+1] = m[i, j] + ds * ms

    return np.concatenate((n[i, N-1] - nL, m[i, N-1] - mL))

def dynamicIVP(G):
    global n, m, p, R, q, w, vs, us, v, u, vt, ut, qt, wt, vst, ust
    n[i] = G[:3]
    m[i] = G[3:]

    for j in range(N-1):
        ps, Rs, ns, ms, qs, ws, vs[i, j], us[i, j], v[i, j], u[i, j], vt[i, j], ut[i, j], qt[i, j], wt[i, j], vst[i, j], ust[i, j] = dynamicODE(p[i, j], R[i, j], n[i, j], m[i, j], q[i, j], w[i, j])
        p[i, j+1] = p[i, j] + ds * ps
        R[i, j+1] = R[i, j] + ds * Rs
        n[i, j+1] = n[i, j] + ds * ns
        m[i, j+1] = m[i, j] + ds * ms
        q[i, j+1] = q[i, j] + ds * qs
        w[i, j+1] = w[i, j] + ds * ws

    return np.concatenate((n[i, N-1] - nLL, m[i, N-1] - mLL))

def staticODE(p, R, n, m):
    v = np.linalg.solve(Kse, R.T @ n + vstar(ds * (j-1)))
    u = np.linalg.solve(Kbt, R.T @ m + ustar(ds * (j-1)))

    ptsb = hat(u) @ rt1 + v
    Tt = 0
    At = -Tt / np.linalg.norm(ptsb)**3 * hat(ptsb) @ hat(ptsb)
    B = hat(rt1) @ At
    Gt = -At @ hat(rt1)
    H = -B @ hat(rt1)
    a = At @ (hat(u) @ ptsb)
    b = hat(rt1) @ a
    d = Kse @ vsstar(ds * (j-1)) - hat(u) @ Kse @ (v - vstar(ds * (j-1))) - a - R.T @ rho * A * g
    c = Kbt @ usstar(ds * (j-1)) - hat(u) @ Kbt @ (u - ustar(ds * (j-1))) - hat(v) @ Kse @ (v - vstar(ds * (j-1))) - b

    Mat = np.block([[Kse + At, Gt], [B, Kbt + H]])
    vs = np.linalg.solve(Mat, np.concatenate((Kbt + H) @ d - Gt @ c, -B @ d + (Kse + At) @ c))

    ps = R @ v
    Rs = R @ hat(u)
    ns = -rho * A * g
    ms = -hat(ps) @ n

    return ps, Rs, ns, ms, vs, us, v, u

def dynamicODE(p, R, n, m, q, w):
    v = np.linalg.solve(Kse + c0 * Bse, R.T @ n + Kse @ vstar(ds * (j-1)) - Bse @ vh[i, j])
    u = np.linalg.solve(Kbt + c0 * Bbt, R.T @ m + Kbt @ ustar(ds * (j-1)) - Bbt @ uh[i, j])
    vt = c0 * v + vh[i, j]
    ut = c0 * u + uh[i, j]
    qt = c0 * q + qh[i, j]
    wt = c0 * w + wh[i, j]

    ptsb1 = hat(u) @ rt1 + v
    At1 = -Tt1 / np.linalg.norm(ptsb1)**3 * hat(ptsb1) @ hat(ptsb1)
    Gt1 = -At1 @ hat(rt1)
    H1 = hat(rt1) @ Gt1
    a1 = At1 @ (hat(u) @ ptsb1)
    b1 = hat(rt1) @ a1

    ptsb2 = hat(u) @ rt2 + v
    At2 = -Tt2 / np.linalg.norm(ptsb2)**3 * hat(ptsb2) @ hat(ptsb2)
    Gt2 = -At2 @ hat(rt2)
    H2 = hat(rt2) @ Gt2
    a2 = At2 @ (hat(u) @ ptsb2)
    b2 = hat(rt2) @ a2

    a = a1 + a2
    b = b1 + b2
    At = At1 + At2
    Gt = Gt1 + Gt2
    H = H1 + H2

    LamdaN = -a + rho * A * (hat(w) @ q + qt) + C @ q * np.abs(q) - R.T @ rho * A * g
    LamdaM = -b + rho * (hat(w) @ J @ w + J @ wt) - hat(v) @ (Kse @ (v - vstar(ds * (j-1))) + Bse @ vt)
    GammaV = hat(u) @ (Kse @ (v - vstar(ds * (j-1))) + Bse @ vt) - Kse @ vsstar(ds * (j-1)) + Bse @ vsh[i, j]
    GammaU = hat(u) @ (Kbt @ (u - ustar(ds * (j-1))) + Bbt @ ut) - Kbt @ usstar(ds * (j-1)) + Bbt @ ush[i, j]

    Mat = np.block([[Kse + c0 * Bse + At, Gt], [Gt.T, Kbt + c0 * Bbt + H]])
    us = np.linalg.solve(Mat, np.concatenate((-Gt.T @ (-GammaV + LamdaN) + (Kse + c0 * Bse + At) @ (-GammaU + LamdaM), (Kbt + c0 * Bbt + H) @ (-GammaV + LamdaN) - Gt @ (-GammaU + LamdaM))))

    vst = c0 * vs + vsh[i, j]
    ust = c0 * us + ush[i, j]

    pts1 = R @ hat(u) @ rt1 + R @ v
    pts2 = R @ hat(u) @ rt2 + R @ v

    nLL = -Tt1 * pts1 / np.linalg.norm(pts1) - Tt2 * pts2 / np.linalg.norm(pts2)
    mLL = -Tt1 * hat(R @ rt1) @ pts1 / np.linalg.norm(pts1) - Tt2 * hat(R @ rt2) @ pts2 / np.linalg.norm(pts2)

    ps = R @ v
    Rs = R @ hat(u)
    ns = rho * A * R @ (hat(w) @ q + qt) - R @ (a + At @ vs + Gt @ us)
    ms = -R @ (b + Gt.T @ vs + H @ us) + rho * R @ (hat(w) @ J @ w + J @ wt) - hat(ps) @ n
    qs = vt - hat(u) @ q + hat(w) @ v
    ws = ut - hat(u) @ w

    return ps, Rs, ns, ms, qs, ws, vs, us, v, u, vt, ut, qt, wt, vst, ust

def visualize1():
    global p, x, y, z, L, N, r, Tt1, Tt2
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)
    for j in range(N):
        x[j] = p[i, j][0]
        y[j] = p[i, j][1]
        z[j] = p[i, j][2]
    filename = f'csvfiles/SingleSectionCR_L{L:.2f}_N{N}_r{r:.4f}_Tt1{Tt1:.2f}_Tt2{Tt2:.2f}.csv'
    print(filename)
    np.savetxt(filename, np.column_stack((x, y, z)), delimiter=',')

def visualize():
    global p, x, y, z, L, N
    if i % 10 == 0:
        x = np.zeros(N)
        y = np.zeros(N)
        z = np.zeros(N)
        for j in range(N):
            x[j] = p[i, j][0]
            y[j] = p[i, j][1]
            z[j] = p[i, j][2]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D(z, y, x)
        ax.set_xlim([-0.05 * L, 1.1 * L])
        ax.set_ylim([-0.1 * L, 0.1 * L])
        ax.set_zlim([-0.05 * L, 0.1 * L])
        ax.set_xlabel('z (m)')
        ax.set_ylabel('y (m)')
        ax.set_zlabel('x (m)')
        plt.grid(True)
        plt.draw()
        plt.pause(0.05)

if __name__ == "__main__":
    SingleSectionCR(Force=10, length=1, resolution=100)