import numpy as np
import astropy.io.fits as pf
from astropy.time import Time
from astropy import units as u
import pandas as pd
# import sys
# sys.path.insert(0, '../ve/vlbi_errors')
# sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
# from components import CGComponent

mas_to_rad = u.mas.to(u.rad)


def kernel(a, b, amp, scale):
    sqdist = np.sum(a**2, 1).reshape(-1, 1) + np.sum(b**2, 1) - 2*np.dot(a, b.T)
    return amp**2*np.exp(-0.5 * (1/(scale*scale)) * sqdist)


def gp_pred(amp, scale, v, t):
    K_ss = kernel(t.reshape(-1, 1), t.reshape(-1, 1), amp, scale)
    L = np.linalg.cholesky(K_ss+1e-6*np.eye(len(t)))
    return np.dot(L, v.reshape(-1, 1))[:, 0]


def gaussian_circ_ft(flux, dx, dy, bmaj, uv):
    """
    FT of circular gaussian at ``uv`` points.

    :param flux:
        Full flux of gaussian.
    :param dx:
        Distance from phase center [mas].
    :param dy:
        Distance from phase center [mas].
    :param bmaj:
        FWHM of a gaussian [mas].
    :param uv:
        2D numpy array of (u,v)-coordinates (dimensionless).
    :return:
        Tuple of real and imaginary visibilities parts.
    """
    shift = [dx*mas_to_rad, dy*mas_to_rad]
    result = np.exp(-2.0*np.pi*1j*(uv @ shift))
    c = (np.pi*bmaj*mas_to_rad)**2/(4. * np.log(2.))
    b = uv[:, 0]**2 + uv[:, 1]**2
    ft = flux*np.exp(-c*b)
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


# TODO: Do not need times_amp & times_phase
def create_data_file(uvfits, outfile, step_amp=60, step_phase=None, use_scans_for_amplitudes=True):
    """
    :param uvfits:
        Path to UVFITS file.
    :param outfile:
        Path to output txt-file.
    :param step_amp: (optional)
        Time interval for constant gain amplitudes [s]. If ``None`` and
        ``use_scans_for_amplitudes`` is ``False`` - use different gain
        amplitude for each time stamp. (default: ``60``)
    :param step_phase:
        Time interval for constant gain phases [s]. If ``None`` - use
        different gain phase for each time stamp. (default: ``None``)
    :param use_scans_for_amplitudes: (optional)
        Boolean. Should we use constant gain amplitude over scan.
        (default: ``True``)
    :return:
    """
    hdus = pf.open(uvfits)
    data = hdus[0].data
    header = hdus[0].header
    freq = header["CRVAL4"]

    df = pd.DataFrame(columns=["times", "ant1", "ant2", "u", "v", "vis_re", "vis_im",
                               "error"])

    for group in data:
        time = Time(group['DATE'] + group['_DATE'], format='jd')
        baseline = group["BASELINE"]
        ant1 = int(baseline//256)
        ant2 = int(baseline-ant1*256)

        try:
            u = group["UU"]
            v = group["VV"]
        except KeyError:
            u = group["UU--"]
            v = group["VV--"]
        if abs(u) < 1.:
            u *= freq
            v *= freq
        data = group["DATA"].squeeze()
        weights = data[:, 2]
        mask = weights <= 0
        if np.alltrue(weights <= 0):
            continue
        weights = np.ma.array(weights, mask=mask)
        weight = np.ma.sum(weights)
        error = 1/np.sqrt(weight)
        vis_re = np.ma.array(data[:, 0], mask=mask)
        vis_im = np.ma.array(data[:, 1], mask=mask)
        vis_re = np.ma.mean(vis_re)
        vis_im = np.ma.mean(vis_im)
        df_ = pd.Series({"times": time, "ant1": ant1, "ant2": ant2, "u": u, "v": v,
                         "vis_re": vis_re, "vis_im": vis_im, "error": error})
        df = df.append(df_, ignore_index=True)

    df["times"] -= np.min(df["times"])
    df["times"] = [dt.sec for dt in df["times"]]

    # All antennas numbers
    antennas = sorted(set(df["ant1"].values.tolist() + df["ant2"].values.tolist()))

    # Time of vis measurements at each antenna
    antennas_times = dict()
    for ant in antennas:
        times = sorted(list(set(df.query("ant1 == @ant | ant2 == @ant")["times"])))
        antennas_times.update({ant: times})
        print("For antenna {} times are: ".format(ant))
        print(times)

    # Indexes of vis measurement at each antenna in unique time stamps
    antennas_times_idx = dict()
    unique_times = sorted(df["times"].unique())
    for ant in antennas:
        antennas_times_idx[ant] = list()
        for t in antennas_times[ant]:
            antennas_times_idx[ant].append(unique_times.index(t))

    # Find scans for each antenna, their central time and indexes of time stamps
    scans_borders = np.where(np.diff(unique_times) > 100)[0]
    print("Data set has {} scans".format(len(scans_borders)+1))
    # Array with scan numbers for each unique time stamp
    scans_idx = np.zeros(len(unique_times), dtype=int)
    scans_idx[: scans_borders[0]+1] = 0
    last_border = scans_borders[0]+1
    i = 1
    for border in scans_borders[1:]:
        scans_idx[last_border:border+1] = i
        i += 1
        last_border = border+1
    scans_idx[last_border:] = i

    # Dictionary with keys - antenna numbers and values - indexes of scans for each time stamp for given antenna
    antennas_scans_idx = dict()
    for ant in antennas:
        antennas_scans_idx[ant] = scans_idx[antennas_times_idx[ant]]

    # Dictionary with keys - scan numbers and values - time of corresponding scan center
    scans = np.unique(scans_idx)
    scans_times = dict()
    for scan in scans:
        scans_times[scan] = np.mean(np.array(unique_times)[scans_idx == scan])
    print("Scan time centers : ", scans_times)

    # Indexes of corresponding visibility measurements in time series for given
    # antenna
    df["id_ant1"] = [antennas_times[ant].index(t) for (ant, t) in
                     df[["ant1", "times"]].values]
    df["id_ant2"] = [antennas_times[ant].index(t) for (ant, t) in
                     df[["ant2", "times"]].values]

    antennas_times_lengths = dict()
    for ant in antennas:
        antennas_times_lengths.update({ant: len(antennas_times[ant])})
        print("For antenna {} number of time measurements is: ".format(ant))
        print(len(antennas_times[ant]))
    # df["ant1_ntimes"] = [antennas_times_lengths[ant] for ant in df["ant1"]]
    # df["ant2_ntimes"] = [antennas_times_lengths[ant] for ant in df["ant2"]]

    # Here using common time grid for re-gridding
    times = df["times"].unique()
    tmin = np.min(times)
    tmax = np.max(times)
    dt = tmax-tmin

    if use_scans_for_amplitudes:
        df["idx_amp_ant1"] = [antennas_scans_idx[ant][idx] for (ant, idx) in
                              df[["ant1", "id_ant1"]].values]
        df["idx_amp_ant2"] = [antennas_scans_idx[ant][idx] for (ant, idx) in
                              df[["ant2", "id_ant2"]].values]
        df["times_amp"] = [scans_times[scan] for scan in df["idx_amp_ant1"].values]

    elif step_amp is not None and not use_scans_for_amplitudes:
        # Re-grid time measurements
        new_idx_dict = {}
        new_times_dict = {}

        # Step size for amplitudes
        step = step_amp
        n = int(dt/step)
        tsamples, step = np.linspace(tmin, tmax, n, retstep=True)
        hist, bins_enges = np.histogram(times, tsamples)
        non_empty = np.where(hist > 0)[0]
        # New times of samples - for amplitudes
        new_times = (tsamples + step/2)[:-1][non_empty]
        new_idx = [a*[i] for i, a in enumerate(hist[non_empty])]
        new_idx = np.array([item for sublist in new_idx for item in sublist], dtype=int)

        # To keep old code working
        for ant in antennas:
            new_idx_dict[ant] = new_idx
            new_times_dict[ant] = new_times
        df["idx_amp_ant1"] = [new_idx_dict[ant][idx] for (ant, idx) in
                              df[["ant1", "id_ant1"]].values]
        df["idx_amp_ant2"] = [new_idx_dict[ant][idx] for (ant, idx) in
                              df[["ant2", "id_ant2"]].values]
        df["times_amp"] = [new_times_dict[ant][idx] for (ant, idx) in
                           df[["ant1", "idx_amp_ant1"]].values]

    else:
        df["idx_amp_ant1"] = df["id_ant1"]
        df["idx_amp_ant2"] = df["id_ant2"]
        df["times_amp"] = df["times"]

    if step_phase is not None:
        # Re-grid time measurements
        new_idx_dict = {}
        new_times_dict = {}

        # Step size for phases
        step = step_phase
        n = int(dt/step)
        tsamples, step = np.linspace(tmin, tmax, n, retstep=True)
        hist, bins_enges = np.histogram(times, tsamples)
        non_empty = np.where(hist > 0)[0]
        # New times of samples - for phases
        new_times = (tsamples+step/2)[:-1][non_empty]
        new_idx = [a*[i] for i, a in enumerate(hist[non_empty])]
        new_idx = np.array([item for sublist in new_idx for item in sublist], dtype=int)

        # To keep old code working
        for ant in antennas:
            new_idx_dict[ant] = new_idx
            new_times_dict[ant] = new_times
        df["idx_phase_ant1"] = [new_idx_dict[ant][idx] for (ant, idx) in
                                df[["ant1", "id_ant1"]].values]
        df["idx_phase_ant2"] = [new_idx_dict[ant][idx] for (ant, idx) in
                                df[["ant2", "id_ant2"]].values]
        df["times_phase"] = [new_times_dict[ant][idx] for (ant, idx) in
                             df[["ant1", "idx_phase_ant1"]].values]
    else:
        df["idx_phase_ant1"] = df["id_ant1"]
        df["idx_phase_ant2"] = df["id_ant2"]
        df["times_phase"] = df["times"]

    df = df[["times",
             "ant1", "ant2",
             "u", "v",
             "vis_re", "vis_im", "error",
             "times_amp", "idx_amp_ant1", "idx_amp_ant2",
             "times_phase", "idx_phase_ant1", "idx_phase_ant2"]]

    df["ant1"] = df["ant1"].astype(int)
    df["ant2"] = df["ant2"].astype(int)
    df.to_csv(outfile, sep=" ", index=False, header=True)
    return df


def add_noise(df, use_global_median_noise=False, use_per_baseline_median_noise=True):
    """
    Add noise as specified in ``error`` columns to ``vis_re`` and ``vis_im`` columns.

    :param df:
        DataFrame with columns ``error``, ``vis_re`` and ``vis_im``.
    :param use_global_median_noise: (optional)
        Noise is estimated using median over all visibilities. Usefull if template uv-data
        has shitty baselines with high noise. (default: ``None``)
    :param use_per_baseline_median_noise: (optional)
        Noise is estimated using median over individual baselines. (default: True)
    :return:
        New DataFrame with noise added.
    """
    df_ = df.copy()

    if use_global_median_noise and use_per_baseline_median_noise:
        raise Exception("Conflicting options of noise estimation!")

    if use_per_baseline_median_noise:
        df_["bl"] = 256*df["ant1"]+df["ant2"]
        g = df_.groupby("bl")
        print("Number of visibilities per-baseline :")
        print(g.bl.count())
        # median error for each baseline
        baselines_med_errors = df_.groupby(["bl"]).median()["error"]
        print("Median error for each baseline :")
        print(baselines_med_errors)
        df_["bl_med_error"] = df_["bl"].map(lambda bl: baselines_med_errors[bl])
        df_["vis_re"] += df_["bl_med_error"].map(lambda x: np.random.normal(0, x, size=1)[0])
        df_["vis_im"] += df_["bl_med_error"].map(lambda x: np.random.normal(0, x, size=1)[0])
        # remove unneeded columns
        del df_["bl"]
        del df_["bl_med_error"]
    elif use_global_median_noise:
        error = df_["error"].median()
        df_["vis_re"] += np.random.normal(0, error, size=df.shape[0])
        df_["vis_im"] += np.random.normal(0, error, size=df.shape[0])
    # Use individual visibilities noise estimated
    else:
        df_["vis_re"] += df_["error"].map(lambda x: np.random.normal(0, x, size=1)[0])
        df_["vis_im"] += df_["error"].map(lambda x: np.random.normal(0, x, size=1)[0])
    return df_


def inject_gains(df_original, outfname=None,
                 amp_gpamp=np.exp(-3), amp_gpphase=np.exp(-3),
                 scale_gpamp=np.exp(6), scale_gpphase=np.exp(5)):
    """
    Inject gains and create new DataFrame with updated visibilities ``df_updated`` and
    dictionary with gains.

    :param df_original:
        Original df with visibilities info.
    :param outfname: (optional)
        Name of file where to save updated DataFrame. If ``None`` than do not save.
        (default: ``None``)
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

    if outfname is not None:
        df_updated.to_csv(outfname, sep=" ", index=False, header=True)

    return df_updated, antennas_gp


def plot_model_gains(gains_dict, savefn=None):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(gains_dict), 2, sharex=True, figsize=(8, 20))
    for i, ant in enumerate(gains_dict):
        amp = gains_dict[ant]["amp"]
        phase = gains_dict[ant]["phase"]
        axes[i, 0].plot(amp.keys(), amp.values(), '.')
        axes[i, 1].plot(phase.keys(), phase.values(), '.')
        axes[i, 1].yaxis.set_ticks_position("right")
    axes[0, 0].set_title("Amplitudes")
    axes[0, 1].set_title("Phases")
    axes[i, 0].set_xlabel("time, s")
    axes[i, 1].set_xlabel("time, s")
    if savefn:
        fig.savefig(savefn, bbox_inches="tight", dpi=300)
    fig.show()
    return fig


# TODO: Finish me
def plot_sampled_gains(posterior_sample, df, fig=None):
    """
    :param posterior_sample:
        DNest file with posterior.
    :param df:
        DataFrame with indexes of gains amplitudes and phases.
    :return:
        Figure.
    """
    with open(posterior_sample, "r") as fo:
        line = fo.readline()

    line = line.split()[1:]
    n_antennas = line.count("ant0")
    gains_len = dict()
    i = 0
    while i < n_antennas:
        gains_len[i] = dict()
        gains_len[i]["amp"] = line.index("phase0") - line.index("amp0")
        line = line[line.index("phase0"):]
        gains_len[i]["phase"] = line.index("amp0") - line.index("phase0")
        line = line[line.index("amp0"):]
        i += 1


def radplot(df, fig=None, color=None, label=None, style="ap"):
    uv = df[["u", "v"]].values
    r = np.hypot(uv[:, 0], uv[:, 1])

    if style == "ap":
        value1 = np.hypot(df["vis_re"].values, df["vis_im"].values)
        value2 = np.arctan2(df["vis_im"].values, df["vis_re"].values)
    elif style == "reim":
        value1 = df["vis_re"].values
        value2 = df["vis_im"].values
    else:
        raise Exception("style can be ap or reim")

    import matplotlib.pyplot as plt
    if fig is None:
        fig, axes = plt.subplots(2, 1, sharex=True)
    else:
        axes = fig.axes

    if color is None:
        color = "#1f77b4"

    axes[0].plot(r, value1, '.', color=color, label=label)
    axes[1].plot(r, value2, '.', color=color)
    if style == "ap":
        axes[0].set_ylabel("Amp, Jy")
        axes[1].set_ylabel("Phase, rad")
    else:
        axes[0].set_ylabel("Re, Jy")
        axes[1].set_ylabel("Im, Jy")
    axes[1].set_xlabel(r"$r_{\rm uv}$, $\lambda$")
    if label is not None:
        axes[0].legend()
    plt.tight_layout()
    fig.show()
    return fig


if __name__ == "__main__":
    # uvfits_fname = "/home/ilya/github/DNest4/code/Examples/UV/J2001+2416_K_2006_06_11_yyk_vis.fits"
    uvfits_fname = "/home/ilya/github/DNest4/code/Examples/UV/0716+714.u.2013_08_20.uvf"
    fname = "/home/ilya/github/bsc/test.txt"

    df = create_data_file(uvfits_fname, fname, step_amp=120, step_phase=30, use_scans_for_amplitudes=True)
# # Load data frame
    # columns = ["times",
    #            "ant1", "ant2",
    #            "u", "v",
    #            "vis_re", "vis_im", "error",
    #            "times_amp", "idx_amp_ant1", "idx_amp_ant2",
    #            "times_phase", "idx_phase_ant1", "idx_phase_ant2"]
    # df = pd.read_csv(fname, sep=" ", names=columns)

    # Zero observed data
    df["vis_re"] = 0
    df["vis_im"] = 0
    # Add model
    re, im = gaussian_circ_ft(flux=3.0, dx=0.0, dy=0.0, bmaj=0.25, uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im

    # Plot
    fig = radplot(df, label="Sky Model")

    # Add gains
    fname_gains = "/home/ilya/github/bsc/test_gains.txt"
    df_updated, gains_dict = inject_gains(df, outfname=fname_gains,
                                          scale_gpamp=np.exp(7),
                                          scale_gpphase=np.exp(5),
                                          amp_gpamp=np.exp(-3),
                                          amp_gpphase=np.exp(-2))
    # Add noise
    df_updated = add_noise(df_updated, use_global_median_noise=True, use_per_baseline_median_noise=False)

    fig = radplot(df_updated, color="#ff7f0e", fig=fig, label="With gains")
    fig = plot_model_gains(gains_dict)

    # import json
    # # Save gains
    # with open("/home/ilya/github/bsc/gains.json", "w") as fo:
    #     json.dump(str(gains_dict), fo)
    # # Load gains
    # with open("/home/ilya/github/bsc/gains_0716.json", "r") as fo:
    #     loaded_gains_dict = json.load(fo)
    # import ast
    # loaded_gains_dict = ast.literal_eval(loaded_gains_dict)
