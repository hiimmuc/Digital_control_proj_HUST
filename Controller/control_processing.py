import matplotlib.pyplot as plt
import numpy as np
import sympy
# from scipy.signal import place_poles as place
from control.matlab import place
from sympy import Symbol
from tqdm import tqdm

MAX_TIME = 20  # seconds
MIN_TIME = 0  # seconds
sampling_rate = 60
Ts = 1/sampling_rate


t = np.linspace(MIN_TIME, MAX_TIME, MAX_TIME*sampling_rate)
l = 0.1  # m length of robot
d = 0.05  # m distance between wheels/2 (center to center)

w_max = np.pi * 2  # rad/s
v_max = 0.8  # m/s
noise = 0  # Gaussian noise

des_poles = np.array([[-2 - 0.5j], [-2 + 0.5j]])

# kinematic transpose matrix
J = np.array([[1, -1, -l+d],
              [1, 1, -l+d],
              [1, -1, l+d],
              [1, 1, l+d]])
J_plus = np.linalg.inv(J.conj().T @ J) @ J.conj().T

# reference trajectory
freq = 2 * np.pi * (1/15)  # rad/s
# t = Symbol('t')
x_ref = np.cos(freq * t)
y_ref = np.sin(freq * t)

# dx_ref = x_ref.diff(t)
# dy_ref = y_ref.diff(t)

# ddx_ref = dx_ref.diff(t)
# ddy_ref = dy_ref.diff(t)
dx_ref = - freq * np.sin(freq * t)
dy_ref = freq * np.cos(freq * t)

ddx_ref = - freq**2 * np.cos(freq * t)
ddy_ref = - freq**2 * np.sin(freq * t)

# circle
q_ref = np.array([x_ref, y_ref, np.arctan2(dy_ref, dx_ref)])
u_ref = np.array([ddx_ref, ddy_ref])

q = np.array([0, 0, -np.pi])
z1 = np.array([q[0], dx_ref[0]], dtype=np.float64).reshape(-1, 1)  # % x, dx state
z2 = np.array([q[1], dy_ref[0]], dtype=np.float64).reshape(-1, 1)  # % y, dy state

v = float(np.sqrt(np.array(z1[1]**2 + z2[1]**2, dtype=np.float64)))  # % body velocity

A = np.array([[0, 1],
              [0, 0]])
B = np.array([[0],
              [1]])
C = np.array([1, 0])


K = np.array(place(A, B, des_poles))

print(f"K matrix is: {K}, shape{K.shape}")

v_r = np.zeros([4, len(t)])  # wheel velocity
q_r = np.zeros([3, len(t)])  # position
w_r = np.zeros([1, len(t)])  # angular velocity
dq_r = np.zeros([3, len(t)])  # body motion (velocity x,y, angular velocity)
D_r = np.zeros([1, len(t)])  # distance from body to destination point


def wrap2pi(p):
    if p == np.pi:
        return ((-p + np.pi) % (2.0 * np.pi) - np.pi) * -1.0
    elif p != np.pi:
        return (p + np.pi) % (2 * np.pi) - np.pi


for k in tqdm(range(0, len(t))):
    z_ref1 = np.array([x_ref[k], dx_ref[k]]).reshape(-1, 1)
    z_ref2 = np.array([y_ref[k], dy_ref[k]]).reshape(-1, 1)
    # print(z_ref1.shape, z_ref2.shape)
    ez1 = z_ref1 - z1
    ez2 = z_ref2 - z2
    # print(ez1.shape, ez2.shape)

    uu = np.array([ddx_ref[k], ddy_ref[k]], dtype='object') + np.array([K @ ez1, K @ ez2], dtype='object').squeeze()
    D = np.linalg.norm([z1[0]-z_ref1[0], z2[0]-z_ref2[0]])

    F = np.array([[np.cos(q[2]), -v * np.sin(q[2])],
                  [np.sin(q[2]), v * np.cos(q[2])]]).squeeze()

    F = F.astype(np.float64)
    uu = uu.reshape(-1, 1).astype(np.float64)

    # print(f"F matrix is \n{F}, {F.shape}\n",
    #       f"uu matrix is\n {uu}, {uu.shape}")

    assert F.shape == (2, 2)
    assert uu.shape == (2, 1)

    vv = np.linalg.solve(F, uu)
    # print('vv res is\n', vv, vv.shape)

    v = v + Ts * float(vv[0])
    u = np.asarray([v, float(vv[1])], dtype='object')
    # print(v, '\n', u, u.shape)

    if abs(u[1]) > w_max:
        u[1] = w_max * np.sign(u[1])
    if abs(v) > v_max:
        v = v_max * np.sign(v)

    dq = np.array([float(u[0] * np.cos(q[2])), float(u[0] * np.sin(q[2])), float(u[1])])

    v_m = J @ dq

    # simulation

    q = q + Ts * dq + np.random.normal(3, 1) * noise
    q[2] = wrap2pi(q[2])

    v_r[:, k] = v_m
    q_r[:, k] = q
    dq_r[:, k] = dq
    D_r[:, k] = D

    z1 = np.array([q[0], u[0]*np.cos(q[2])]).reshape(-1, 1)
    z2 = np.array([q[1], u[0]*np.sin(q[2])]).reshape(-1, 1)


plt.figure()
plt.plot(q_r[0, :], q_r[1, :])
plt.plot(x_ref, y_ref, 'r--')
plt.title("xy plane charting")
plt.legend(['robot', 'reference'])

plt.figure()
plt.plot(t, v_r[0, :], 'r')
plt.plot(t, v_r[1, :], 'g')
plt.plot(t, v_r[2, :], 'b')
plt.plot(t, v_r[3, :], 'c')

plt.title("Individual wheel input velocity")
plt.legend(['wheel 1', 'wheel 2', 'wheel 3', 'wheel 4'])

plt.figure()
plt.plot(t, D_r.squeeze())
plt.grid(True)
plt.title("Error distance from body to reference point")
plt.show()
plt.close()
