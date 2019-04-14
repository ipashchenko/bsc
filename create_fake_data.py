import numpy as np
import pandas as pd


def kernel(a, b, amp, scale):
    sqdist = np.sum(a**2, 1).reshape(-1, 1) + np.sum(b**2, 1) - 2*np.dot(a, b.T)
    return amp**2*np.exp(-0.5 * (1/(scale*scale)) * sqdist)


def gp_pred(amp, scale, v, t):
    K_ss = kernel(t.reshape(-1, 1), t.reshape(-1, 1), amp, scale)
    L = np.linalg.cholesky(K_ss+1e-6*np.eye(len(t)))
    return np.dot(L, v.reshape(-1, 1))[:, 0]


def inject_gains(df_original,
                 amp_gpamp=np.exp(-3), amp_gpphase=np.exp(-3),
                 scale_gpamp=np.exp(6), scale_gpphase=np.exp(5)):
    """
    Inject gains and create new DataFrame with updated visibilities ``df_updated`` and
    dictionary with gains.

    :param df_original:
        Original df with visibilities info.
    :param amp_gpamp: (optional)
        Amplitude of the GP that describes time dependence of the gain amplitudes.
        (default: ``np.exp(-3)``.

    :return:
        Updated DataFrame and dictionary with gains for each antenna.
    """
    df = df_original.copy()

    # All antennas numbers
    antennas = sorted(set(df["ant1"].values.tolist() + df["ant2"].values.tolist()))

    # Time of vis measurements at each antenna
    antennas_times = dict()
    for ant in antennas:
        times = sorted(list(set(df.query("ant1 == @ant | ant2 == @ant")["times"])))
        antennas_times.update({ant: times})

    # For each antenna - create GP of amp and phase
    antennas_gp = dict()
    for ant in antennas:
        ant_time = antennas_times[ant]
        # Amplitude
        v = np.random.normal(0, 1, size=len(ant_time))
        amp = 1 + gp_pred(amp_gpamp, scale_gpamp, v, np.array(ant_time))
        # Phase
        v = np.random.normal(0, 1, size=len(ant_time))
        phase = gp_pred(amp_gpphase, scale_gpphase, v, np.array(ant_time))
        antennas_gp[ant] = {"amp": {t: a for (t, a) in zip(ant_time, amp)},
                            "phase": {t: p for (t, p) in zip(ant_time, phase)}}

    # Now calculate visibilities
    vis_real_gained = list()
    vis_imag_gained = list()
    for t, ant1, ant2, vis_real, vis_imag in df[['times', 'ant1', 'ant2', 'vis_re', 'vis_im']].values:
        amp1 = antennas_gp[ant1]["amp"][t]
        amp2 = antennas_gp[ant2]["amp"][t]
        phase1 = antennas_gp[ant1]["phase"][t]
        phase2 = antennas_gp[ant2]["phase"][t]
        vis_real_new = amp1*amp2*(np.cos(phase1-phase2)*vis_real-np.sin(phase1-phase2)*vis_imag)
        vis_imag_new = amp1*amp2*(np.cos(phase1-phase2)*vis_imag+np.sin(phase1-phase2)*vis_real)
        vis_real_gained.append(vis_real_new)
        vis_imag_gained.append(vis_imag_new)

    # Create updated df that can be used for inference
    df_updated = df_original.copy()
    df_updated["vis_re"] = vis_real_gained
    df_updated["vis_im"] = vis_imag_gained

    return df_updated, antennas_gp


if __name__ == "__main__":
    # Load data frame
    columns = ["times",
               "ant1", "ant2",
               "u", "v",
               "vis_re", "vis_im", "error",
               "times_amp", "idx_amp_ant1", "idx_amp_ant2",
               "times_phase", "idx_phase_ant1", "idx_phase_ant2"]
    fname = "/home/ilya/CLionProjects/ve2/uv_data.txt"
    df = pd.read_csv(fname, sep=" ", names=columns)

    df_updated, gains_dict = inject_gains(df)