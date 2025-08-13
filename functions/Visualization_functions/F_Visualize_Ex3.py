import numpy as np
import matplotlib.pyplot as plt

def F_Visualize_Ex3(Y_ens, S_obs, S_true, S_ens, xPlus_3d, Sll_2d, Cxx_3d, Cxpxp_3d):
    clr_true = 'k'
    marker_true = '.'
    line_true = 'k-'
    msize_true = 15

    clr_obs = 'k'
    marker_obs = '^'
    line_obs = 'k-'
    msize_obs = 4

    clr_mod = '#d10000'
    marker_mod = 'o'
    line_mod = '--'
    clr_mod2 = '#d10000'

    clr_upd = '#527dea'
    marker_upd = 'o'
    line_upd = ':'
    clr_upd2 = 'w'

    clr_ens = (0.8, 0.8, 0.8)
    line_ens = '-'

    lw = 2
    msize = 6
    fsize = 12
    opt_save = True
    opt_saveData = False

    Ne = Y_ens.shape[1]
    sim = [1, 24]  # simulation phase of model

    # ---- Visualizations ----

    # truth and observations
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    for nn in range(Ne):
        plt.plot(Y_ens[:, nn], line_ens, color=clr_ens, linewidth=lw)
    h2, = plt.plot(S_obs, line_obs, marker=marker_obs, markeredgecolor=clr_obs,
                   markerfacecolor=clr_obs, linewidth=lw, markersize=msize)
    h3, = plt.plot(S_true, line_true, marker=marker_true, markeredgecolor=clr_true,
                   markerfacecolor=clr_true, linewidth=lw, markersize=msize_true)
    plt.legend([h3, h2], ['Truth', 'Observation'])
    plt.xlabel('time')
    plt.xlim([-0.5, (sim[1] - sim[0] + 1) + 0.5])
    plt.ylabel('S')
    plt.title('Ground truth and observation ensemble')
    plt.grid(True)
    if opt_save:
        plt.savefig('figures/Ex3_Obs.jpg', dpi=150)
    plt.show()

    # model simulation (open loop)
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    for nn in range(Ne):
        plt.plot(S_ens[:, nn], line_ens, color=clr_ens, linewidth=lw)
    h2, = plt.plot(np.mean(S_ens, axis=1), line_mod, marker=marker_mod, markeredgecolor=clr_mod,
                   markerfacecolor=clr_mod2, linewidth=lw, color=clr_mod, markersize=msize)
    plt.legend([h2], ['Model mean'])
    plt.xlabel('time')
    plt.xlim([-0.5, (sim[1] - sim[0] + 1) + 0.5])
    plt.ylabel('S')
    plt.title('Water storage estimates - OL')
    plt.grid(True)
    if opt_save:
        plt.savefig('figures/Ex3_OL.jpg', dpi=150)
    plt.show()

    # DA results
    length_sim = sim[1] - sim[0] + 1
    xPlus_ensMean = np.mean(xPlus_3d[0, :, :], axis=0)

    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    for nn in range(Ne):
        plt.plot(xPlus_3d[0, nn, :], line_ens, color=clr_ens, linewidth=lw)
    h2, = plt.plot(xPlus_ensMean, line_upd, marker=marker_upd, markeredgecolor=clr_upd,
                   markerfacecolor=clr_upd, linewidth=lw, color=clr_upd, markersize=msize)
    plt.legend([h2], ['Model mean'])
    plt.xlabel('time')
    plt.xlim([-0.5, length_sim + 0.5])
    plt.ylabel('S')
    plt.title('Water storage estimates - DA')
    plt.grid(True)
    if opt_save:
        plt.savefig('figures/Ex3_DA.jpg', dpi=150)
    plt.show()

    # ensemble prediction
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    plt.plot(np.mean(S_ens, axis=1), line_mod, marker=marker_mod, markeredgecolor=clr_mod,
             markerfacecolor=clr_mod2, linewidth=lw, color=clr_mod, markersize=msize)
    plt.plot(xPlus_ensMean, line_upd, marker=marker_upd, markeredgecolor=clr_upd,
             markerfacecolor=clr_upd2, linewidth=lw, color=clr_upd, markersize=msize)
    plt.plot(S_obs, line_obs, marker=marker_obs, markeredgecolor=clr_obs,
             markerfacecolor=clr_obs, linewidth=lw, markersize=msize)
    plt.plot(S_true, line_true, marker=marker_true, markeredgecolor=clr_true,
             markerfacecolor=clr_true, linewidth=lw, markersize=msize_true)
    plt.legend(['OL', 'DA', 'Observation', 'Truth'])
    plt.xlabel('time')
    plt.xlim([-0.5, length_sim + 0.5])
    plt.ylabel('S')
    plt.title('Water storage estimates (ensemble average)')
    plt.grid(True)
    if opt_save:
        plt.savefig('figures/Ex3_EnsAv.jpg', dpi=150)
    plt.show()

    # empirical standard deviations
    plt.figure()
    plt.gca().tick_params(labelsize=fsize)
    plt.plot(np.sqrt(Sll_2d), line_obs, marker=marker_obs, markeredgecolor=clr_obs,
             markerfacecolor=clr_obs, linewidth=lw, markersize=msize)

    test1 = np.zeros(length_sim)
    test1[:] = Cxx_3d[0, 0, :]
    plt.plot(np.sqrt(test1), line_mod, marker=marker_mod, markeredgecolor=clr_mod,
             markerfacecolor=clr_mod2, linewidth=lw, color=clr_mod, markersize=msize)

    test2 = np.zeros(length_sim)
    test2[:] = Cxpxp_3d[0, 0, :]
    plt.plot(np.sqrt(test2), line_upd, marker=marker_upd, markeredgecolor=clr_upd,
             markerfacecolor=clr_upd2, linewidth=lw, color=clr_upd, markersize=msize)

    plt.xlabel('time')
    plt.xlim([-0.5, length_sim + 0.5])
    plt.ylabel('S sigma')
    plt.title('Uncertainty (based on ensemble spread)')
    plt.grid(True)
    plt.legend(['Observations', 'Model Prediction', 'Model Update'])
    if opt_save:
        plt.savefig('figures/Ex3_Uncert.jpg', dpi=150)
    plt.show()
