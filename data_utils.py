import numpy as np
import astropy.io.fits as pf
from astropy.time import Time
import pandas as pd


def kernel(a, b, amp, scale):
    sqdist = np.sum(a**2, 1).reshape(-1, 1) + np.sum(b**2, 1) - 2*np.dot(a, b.T)
    return amp**2*np.exp(-0.5 * (1/(scale*scale)) * sqdist)


def gp_pred(amp, scale, v, t):
    K_ss = kernel(t.reshape(-1, 1), t.reshape(-1, 1), amp, scale)
    L = np.linalg.cholesky(K_ss+1e-6*np.eye(len(t)))
    return np.dot(L, v.reshape(-1, 1))[:, 0]


def create_data_file(uvfits, outfile, step_amp=60, step_phase=None):
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

        u = group["UU"]
        v = group["VV"]
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

    # Indexes of corresponding visibility measurements in time series for given
    # antenna
    df["id_ant1"] = [antennas_times[ant].index(t) for (ant, t) in
                     df[["ant1", "times"]].values]
    df["id_ant2"] = [antennas_times[ant].index(t) for (ant, t) in
                     df[["ant2", "times"]].values]

    # antennas_times_lengths = dict()
    # for ant in antennas:
    #     antennas_times_lengths.update({ant: len(antennas_times[ant])})
    # df["ant1_ntimes"] = [antennas_times_lengths[ant] for ant in df["ant1"]]
    # df["ant2_ntimes"] = [antennas_times_lengths[ant] for ant in df["ant2"]]


    # Re-grid time measurements
    new_idx_dict = {}
    new_times_dict = {}
    for ant in antennas:
        # Step size for amplitudes
        step = step_amp
        times = np.array(antennas_times[ant])
        tmin = np.min(times)
        tmax = np.max(times)
        dt = tmax-tmin
        n = int(dt/step)
        tsamples, step = np.linspace(tmin, tmax, n, retstep=True)
        hist, bins_enges = np.histogram(times, tsamples)
        non_empty = np.where(hist > 0)[0]
        # New times of samples
        new_times = (tsamples + step/2)[:-1][non_empty]
        new_idx = [a*[i] for i, a in enumerate(hist[non_empty])]
        new_idx = np.array([item for sublist in new_idx for item in sublist], dtype=int)
        new_idx_dict[ant] = new_idx
        new_times_dict[ant] = new_times
    df["idx_amp_ant1"] = [new_idx_dict[ant][idx]  for (ant, idx) in
                          df[["ant1", "id_ant1"]].values]
    df["idx_amp_ant2"] = [new_idx_dict[ant][idx]  for (ant, idx) in
                          df[["ant2", "id_ant2"]].values]
    df["times_amp"] = [new_times_dict[ant][idx] for (ant, idx) in
                       df[["ant1", "idx_amp_ant1"]].values]


    if step_phase is not None:
        # Re-grid time measurements
        new_idx_dict = {}
        new_times_dict = {}
        for ant in antennas:
            # Step size for amplitudes
            step = step_phase
            times = np.array(antennas_times[ant])
            tmin = np.min(times)
            tmax = np.max(times)
            dt = tmax-tmin
            n = int(dt/step)
            tsamples, step = np.linspace(tmin, tmax, n, retstep=True)
            hist, bins_enges = np.histogram(times, tsamples)
            non_empty = np.where(hist > 0)[0]
            # New times of samples
            new_times = (tsamples+step/2)[:-1][non_empty]
            new_idx = [a*[i] for i, a in enumerate(hist[non_empty])]
            new_idx = np.array(
                [item for sublist in new_idx for item in sublist], dtype=int)
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

    df.to_csv(outfile, sep=" ", index=False, header=False)
    return df


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
        df_updated.to_csv(outfname, sep=" ", index=False, header=False)

    return df_updated, antennas_gp


if __name__ == "__main__":
    uvfits_fname = "/home/ilya/github/DNest4/code/Examples/UV/J2001+2416_K_2006_06_11_yyk_vis.fits"
    fname = "/home/ilya/github/bsc/uv_data.txt"

    df = create_data_file(uvfits_fname, fname, step_amp=120)
    # # Load data frame
    # columns = ["times",
    #            "ant1", "ant2",
    #            "u", "v",
    #            "vis_re", "vis_im", "error",
    #            "times_amp", "idx_amp_ant1", "idx_amp_ant2",
    #            "times_phase", "idx_phase_ant1", "idx_phase_ant2"]
    # df = pd.read_csv(fname, sep=" ", names=columns)
    fname_gains = "/home/ilya/github/bsc/uv_data_gains.txt"
    df_updated, gains_dict = inject_gains(df, outfname=fname_gains)

    import json
    # Save gains
    with open("/home/ilya/github/bsc/gains.json", "w") as fo:
        json.dump(str(gains_dict), fo)
    # Load gains
    with open("/home/ilya/github/bsc/gains.json", "r") as fo:
        loaded_gains_dict = json.load(fo)
    import ast
    loaded_gains_dict = ast.literal_eval(loaded_gains_dict)
