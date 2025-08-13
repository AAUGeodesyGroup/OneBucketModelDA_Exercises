import numpy as np
import matplotlib.pyplot as plt
import os

def F_Visualize_Ex1(P, S, S0, Q, K):
    lw = 2
    msize = 15
    fsize = 20
    opt_save = 1

    # Ensure figures directory exists
    if opt_save == 1 and not os.path.exists("figures"):
        os.makedirs("figures")

    # precipitation
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    plt.bar(np.arange(1, len(P) + 1), P, color='k')
    plt.xlabel('time')
    plt.xlim([-0.5, len(P) + 0.5])
    plt.xticks(np.arange(0, 51, 10))
    plt.ylabel('P-E')
    plt.title('Input data: precipitation - evapotranspiration')
    plt.box(True)
    if opt_save == 1:
        plt.savefig('figures/Ex1_P.jpeg', dpi=300)
    plt.show()
    plt.close()

    # storage
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    plt.plot(np.arange(0, len(P) + 1), np.insert(S, 0, S0), '.-k', linewidth=lw, markersize=msize)
    plt.xlabel('time')
    plt.xlim([-0.5, len(P) + 0.5])
    plt.xticks(np.arange(0, 51, 10))
    plt.ylabel('S')
    plt.ylim([0, np.max(np.insert(S, 0, S0)) + 2])
    plt.title('Water storage estimates - single model simulation')
    plt.box(True)
    if opt_save == 1:
        plt.savefig('figures/Ex1_S.jpeg', dpi=300)
    plt.show()
    plt.close()

    # runoff/discharge
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    plt.plot(np.arange(1, len(P) + 1), Q, '.-k', linewidth=lw, markersize=msize)
    plt.xlabel('time')
    plt.xlim([-0.5, len(P) + 0.5])
    plt.xticks(np.arange(0, 51, 10))
    plt.ylabel('R')
    plt.ylim([0, np.max(np.insert(S, 0, S0)) + 2])
    plt.title('River discharge estimates - single model simulation')
    plt.box(True)
    if opt_save == 1:
        plt.savefig('figures/Ex1_R.jpeg', dpi=300)
    plt.show()
    plt.close()

    # relation between S and Q
    S_sim = np.arange(0, 25.1, 0.1)
    Q_sim = K * S_sim
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    plt.plot(S_sim, Q_sim, '-k', linewidth=lw, markersize=msize)
    plt.xlabel('S')
    plt.ylabel('R')
    plt.title('Relation between storage and discharge - single model simulation')
    plt.box(True)
    plt.axis('equal')
    plt.xlim([0, 25])
    plt.ylim([0, 8])
    if opt_save == 1:
        plt.savefig('figures/Ex1_StoR.jpeg', dpi=300)
    plt.show()
    plt.close()
