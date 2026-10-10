import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

                     
def second_derivative_matrix(N, h):
    diagonal = np.eye(N) * -2
    upper_off_diagonal = np.eye(N, k=1)
    lower_off_diagonal = np.eye(N, k=-1)
    D2 = (1/h**2) * (diagonal + upper_off_diagonal + lower_off_diagonal)
    return D2

def hamiltonian(N, h, V):
    D2 = second_derivative_matrix(N, h)
    H = -0.5 * D2
    V_matrix = np.diag(V)
    return H + V_matrix


x_min = -100
x_max = 100
N = 1000
h = (x_max - x_min)/(N + 1)
x = np.linspace(x_min + h, x_max - h, N)

x0 = -50
a = 5.0
k0 = 3.0

psi0 = np.exp(-(x - x0) ** 2 / a ** 2) * np.exp(1j *k0 *x)
norm = np.sqrt(np.sum(np.abs(psi0)**2)*h)
psi = psi0 / norm

if __name__ == "__main__":

    print(np.sum(np.abs(psi)**2)*h)

    plt.plot(x, np.abs(psi)**2, label="|Psi| ^ 2")
    plt.plot(x, np.real(psi), label="Re(Psi)")
    plt.legend()
    plt.show()

    dt = 0.005
    V = np.zeros(N)

    t_final = 1
    n_steps = int(t_final / dt)

    H = hamiltonian(N, h, V)
    A = np.eye(N, dtype=complex) + 0.5j * dt * H
    B = np.eye(N, dtype=complex) - 0.5j * dt * H

    psi_t = psi.copy()
    for step in range(n_steps):
        b = B @ psi_t
        psi_t = np.linalg.solve(A, b)

    print(np.sum(np.abs(psi_t)**2)*h)

    plt.plot(x, np.abs(psi)**2, label="t = 0")
    plt.plot(x, np.abs(psi_t) ** 2, label="t = t final")
    plt.legend()
    plt.title("probability density")
    plt.show()

    x = np.linspace(x_min + h, x_max - h, N)
    V = []
    for xn in x:
        if np.abs(xn)<0.5:
            V.append(20)
        else:
            V.append(0)
    V = np.array(V)
    H = hamiltonian(N, h, V)
    A = np.eye(N, dtype=complex) + 0.5j * dt * H
    B = np.eye(N, dtype=complex) - 0.5j * dt * H

    psi_t = psi.copy()
    for step in range(n_steps):
        b = B @ psi_t
        psi_t = np.linalg.solve(A, b)

    plt.plot(x, np.abs(psi) ** 2, label="t = 0")
    plt.plot(x, np.abs(psi_t) ** 2, label="t = t final")
    plt.legend()
    plt.title("probability density with potential wall")
    plt.show()

    plt.plot(x, np.abs(psi_t) ** 2)
    plt.yscale('log')
    plt.ylim(1e-14, 1)
    plt.axvline(-0.5, color ='k', linestyle='--')
    plt.axvline(0.5, color = 'k', linestyle='--')
    plt.title("log scale - transmitted tail")
    plt.show()


    t_final = 40
    n_steps = int(t_final / dt)

    H = hamiltonian(N, h, V)
    A = np.eye(N, dtype=complex) + 0.5j * dt * H
    B = np.eye(N, dtype=complex) - 0.5j * dt * H

    psi_t = psi.copy()
    for step in range(n_steps):
        b = B @ psi_t
        psi_t = np.linalg.solve(A, b)

    R = np.sum((np.abs(psi_t[x < -0.5]) ** 2)*h )
    T = np.sum((np.abs(psi_t[x > 0.5]) ** 2) * h)
    P_in = np.sum((np.abs(psi_t[np.abs(x) <= 0.5]) ** 2) *h)
    print(R, T, P_in, R + T + P_in)


    t_video = 6
    fps = 25
    n_frames = t_video * fps
    step_width = int(n_steps / n_frames)

    snapshots = []
    psi_t = psi.copy()
    for step in range(n_steps):
        b = B @ psi_t
        psi_t = np.linalg.solve(A, b)
        if step % step_width == 0:
            snapshots.append(np.abs(psi_t) ** 2)

    print(len(snapshots))

    fig, ax = plt.subplots()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 0.2)
    ax.axvline(-0.5, color='k', linestyle='--')
    ax.axvline(0.5, color='k', linestyle='--')
    line, = ax.plot([], [])

    def update(frame):
        line.set_data(x, snapshots[frame])
        return line,

    anim = FuncAnimation(fig, update, frames=len(snapshots), blit=True)
    anim.save('tunnelling.gif', writer=PillowWriter(fps=fps))

