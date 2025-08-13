import numpy as np
import matplotlib.pyplot as plt
import os

def F_Visualize_Ex2(P_ens, S_ens, Q_ens, sim):
    """
    Visualization for Exercise 2 ensemble runs.

    Parameters
    ----------
    P_ens : ndarray
        Precipitation ensemble, shape (T, Ne)
    S_ens : ndarray
        Storage ensemble, shape (T, Ne)
    Q_ens : ndarray
        Discharge ensemble, shape (T, Ne)
    sim : list or tuple
        Simulation indices [start_idx, end_idx] (Python 0-based)
    """

    clr_ens = [0.5, 0.5, 0.5]  # ensemble color
    lw = 2
    msize = 15
    fsize = 12
    opt_save = 1
    Ne = P_ens.shape[1]

    # Ensure figures directory exists
    os.makedirs('figures', exist_ok=True)

    time_range = np.arange(sim[0], sim[1])

    # === Precipitation ===
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    for ii in range(Ne):
        h1, = plt.plot(time_range, P_ens[time_range, ii], '.-', color=clr_ens,
                       linewidth=lw, markersize=msize)
    h2, = plt.plot(time_range, np.mean(P_ens[time_range, :], axis=1), '.-', color='k',
                   linewidth=lw, markersize=msize)
    plt.xlabel('time')
    plt.xlim([-0.5, (sim[1] - sim[0]) + 0.5])
    plt.ylabel('P-E')
    plt.title('Perturbed input data: precipitation - evapotranspiration ensemble')
    plt.legend([h2, h1], ['Ensemble mean', 'Ensemble members'])
    plt.box(True)
    if opt_save == 1:
        plt.savefig('figures/Ex2_P_OL.jpg', dpi=300)

    # === Storage ===
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    for ii in range(Ne):
        h1, = plt.plot(time_range, S_ens[:, ii], '.-', color=clr_ens,
                       linewidth=lw, markersize=msize)
    h2, = plt.plot(time_range, np.mean(S_ens, axis=1), '.-', color='k',
                   linewidth=lw, markersize=msize)
    plt.xlabel('time')
    plt.xlim([-0.5, (sim[1] - sim[0]) + 0.5])
    plt.ylabel('S')
    plt.title('Water storage estimates - model simulation ensemble')
    plt.legend([h2, h1], ['Ensemble mean', 'Ensemble members'])
    plt.ylim([0, np.max(S_ens) + 2])
    plt.box(True)
    if opt_save == 1:
        plt.savefig('figures/Ex2_S_OL.jpg', dpi=300)

    # === Discharge ===
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    for ii in range(Ne):
        h1, = plt.plot(time_range, Q_ens[:, ii], '.-', color=clr_ens,
                       linewidth=lw, markersize=msize)
    h2, = plt.plot(time_range, np.mean(Q_ens, axis=1), '.-', color='k',
                   linewidth=lw, markersize=msize)
    plt.xlabel('time')
    plt.xlim([-0.5, (sim[1] - sim[0]) + 0.5])
    plt.ylabel('R')
    plt.title('River discharge estimates - model simulation ensemble')
    plt.ylim([0, np.max(S_ens) + 2])  # same ylim as MATLAB
    plt.legend([h2, h1], ['Ensemble mean', 'Ensemble members'])
    plt.box(True)
    if opt_save == 1:
        plt.savefig('figures/Ex2_R_OL.jpg', dpi=300)

    plt.show()
