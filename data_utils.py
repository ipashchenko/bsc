import numpy as np
import astropy.io.fits as pf
from astropy.time import Time
from astropy import units as u
from tqdm import tqdm
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from uv_data import UVData

label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size
matplotlib.rcParams['font.size'] = label_size
matplotlib.rcParams['legend.fontsize'] = label_size
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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


# TODO: Check STOKES/IF
def create_data_file(uvfits, outfile, STOKES, IF, step_amp=30, step_phase=30, use_scans_for_amplitudes=False,
                     calculate_noise=False, antennas_to_skip=None):
    """
    :param uvfits:
        Path to UVFITS file.
    :param outfile:
        Path to output txt-file.
    :param STOKES:
        Stokes number (RR = 0, etc.)
    :param IF:
        IF number.
    :param step_amp: (optional)
        Time interval for constant gain amplitudes [s]. If ``None`` and
        ``use_scans_for_amplitudes`` is ``False`` - use different gain
        amplitude for each time stamp. (default: ``30``)
    :param step_phase:
        Time interval for constant gain phases [s]. If ``None`` - use
        different gain phase for each time stamp. (default: ``30``)
    :param use_scans_for_amplitudes: (optional)
        Boolean. Should we use constant gain amplitude over scan.
        (default: ``False``)
    :return:
    """
    hdus = pf.open(uvfits)
    data_all = hdus[0].data
    header = hdus[0].header
    freq = header["CRVAL4"]

    uvdata = UVData(uvfits)
    noise = uvdata.noise(use_V=False)

    df = pd.DataFrame(columns=["times", "ant1", "ant2", "u", "v", "vis_re", "vis_im",
                               "error"])

    suffix = None
    for possible in ['UU', 'UU--', 'UU---SIN']:
        if possible in data_all.parnames:
            suffix = possible
            break
    if suffix is None:
        raise Exception("CHECK BASELINE COORDINATES SYSTEM!")

    for group in tqdm(data_all):
        time = Time(group['DATE'] + group['_DATE'], format='jd')
        baseline = group["BASELINE"]
        ant1 = int(baseline//256)
        ant2 = int(baseline-ant1*256)
        if antennas_to_skip is not None:
            if ant1 in antennas_to_skip or ant2 in antennas_to_skip:
                print("Skipping baseline {} - {}".format(ant1, ant2))
                continue
        bl_noise = noise[float(baseline)]/np.sqrt(2.)

        u = group[suffix]
        v = group[suffix]
        if abs(u) < 1.:
            u *= freq
            v *= freq
        # data = group["DATA"].squeeze()
        # IF, STOKES, COMPLEX
        data = group["DATA"][0, 0, IF, 0, STOKES, :]
        weights = data[..., 2]
        mask = weights <= 0
        if np.alltrue(weights <= 0):
            continue
        weights = np.ma.array(weights, mask=mask)
        weight = np.ma.sum(weights)
        if calculate_noise:
            error = bl_noise[IF, STOKES]
        else:
            error = 1/np.sqrt(weight)
        vis_re = np.ma.array(data[..., 0], mask=mask)
        vis_im = np.ma.array(data[..., 1], mask=mask)
        vis_re = np.ma.mean(vis_re)
        vis_im = np.ma.mean(vis_im)
        df_ = pd.Series({"times": time, "ant1": ant1, "ant2": ant2, "u": u, "v": v,
                         "vis_re": vis_re, "vis_im": vis_im, "error": error})
        df = df.append(df_, ignore_index=True)

    mint = np.min(df["times"])
    df["times"] = df["times"].apply(lambda x: x - mint)
    df["times"] = df['times'].apply(lambda x: x.sec)
    # df["times"] = [dt.sec for dt in df["times"]]

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


def create_data_file_many_IFs(uvfits, outfile, step_amp=30, step_phase=30, use_scans_for_amplitudes=False,
                              calculate_noise=False, antennas_to_skip=None):
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
    data_all = hdus[0].data
    header = hdus[0].header
    freq = header["CRVAL4"]

    uvdata = UVData(uvfits)
    noise = uvdata.noise(use_V=False)

    df = pd.DataFrame(columns=["times", "ant1", "ant2", "u", "v", "STOKES", "IF", "vis_re", "vis_im",
                               "error"])

    suffix = None
    for possible in ['UU', 'UU--', 'UU---SIN']:
        if possible in data_all.parnames:
            suffix = possible
            break
    if suffix is None:
        raise Exception("CHECK BASELINE COORDINATES SYSTEM!")

    for group in tqdm(data_all):
        time = Time(group['DATE'] + group['_DATE'], format='jd')
        baseline = group["BASELINE"]
        ant1 = int(baseline//256)
        ant2 = int(baseline-ant1*256)
        if antennas_to_skip is not None:
            if ant1 in antennas_to_skip or ant2 in antennas_to_skip:
                print("Skipping baseline {} - {}".format(ant1, ant2))
                continue

        bl_noise = noise[float(baseline)]/np.sqrt(2.)

        u = group[suffix]
        v = group[suffix]

        if abs(u) < 1.:
            u *= freq
            v *= freq

        # IF, STOKES, COMPLEX
        data = group["DATA"][0, 0, :, 0, :, :]
        n_IF = data.shape[0]
        n_stokes = data.shape[1]
        for i_IF in range(n_IF):
            for i_stokes in range(n_stokes):
                weight = data[i_IF, i_stokes][2]
                if weight <= 0:
                    continue

                if calculate_noise:
                    error = bl_noise[IF, STOKES]
                else:
                    error = 1/np.sqrt(weight)

                vis_re = data[i_IF, i_stokes][0]
                vis_im = data[i_IF, i_stokes][1]

                df_ = pd.Series({"times": time, "ant1": ant1, "ant2": ant2, "u": u, "v": v, "STOKES": i_stokes, "IF": i_IF,
                                 "vis_re": vis_re, "vis_im": vis_im, "error": error})
                df = df.append(df_, ignore_index=True)

    mint = np.min(df["times"])
    df["times"] = df["times"].apply(lambda x: x - mint)
    df["times"] = df['times'].apply(lambda x: x.sec)

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
             "u", "v", "STOKES", "IF",
             "vis_re", "vis_im", "error",
             "times_amp", "idx_amp_ant1", "idx_amp_ant2",
             "times_phase", "idx_phase_ant1", "idx_phase_ant2"]]

    df["ant1"] = df["ant1"].astype(int)
    df["ant2"] = df["ant2"].astype(int)
    df["IF"] = df["IF"].astype(int)
    df["STOKES"] = df["STOKES"].astype(int)
    df.to_csv(outfile, sep=" ", index=False, header=True)
    return df


def add_noise(df, use_global_median_noise=False, use_per_baseline_median_noise=True, outfname=None):
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

    if outfname is not None:
        df_.to_csv(outfname, sep=" ", index=False, header=True)

    return df_


def inject_gains(df_original, outfname=None, non_zero_mean_phase=False,
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
        if non_zero_mean_phase:
            mean_phase = np.random.uniform(-np.pi, np.pi)
        else:
            mean_phase = 0
        v = np.random.normal(0, 1, size=len(ant_time))
        phase = gp_pred(amp_gpphase, scale_gpphase, v, np.array(ant_time)) + mean_phase
        antennas_gp[ant] = {"amp": {t: a for (t, a) in zip(ant_time, amp)},
                            "mean_phase": mean_phase,
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
    fig, axes = plt.subplots(len(gains_dict), 2, sharex=True, figsize=(8, 20))
    for i, ant in enumerate(gains_dict):
        amp = gains_dict[ant]["amp"]
        phase = gains_dict[ant]["phase"]
        axes[i, 0].plot(list(amp.keys()), list(amp.values()), '.')
        axes[i, 1].plot(list(phase.keys()), list(phase.values()), '.')
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
def process_sampled_gains(posterior_sample, df_fitted, jitter_first=True, n_comp=2, plotfn=None,
                          add_mean_phase=False):
    """
    :param posterior_sample:
        DNest file with posterior.
    :param df_fitted:
        Dataframe fitted (created by ``create_data_file`` and ``inject_gains``).
    :return:
        Dictionary with keys - antenna numbers (as in ``gains_dict``), times, ``amp``, ``phase``
        and values - samples of the posterior distribution for amplitude and phase at given time
        for given antenna.
    """
    samples = np.loadtxt(posterior_sample, skiprows=1)
    first_gain_index = n_comp*4
    if jitter_first:
        first_gain_index += 1
    gains_len = dict()
    gains_post = dict()
    j = first_gain_index
    antennas = set(list(df_fitted.ant1.unique()) + list(df_fitted.ant2.unique()))
    for ant in antennas:
        gains_len[ant] = dict()
        gains_len[ant]["amp"] = len(set(df_fitted.query("ant1 == @ant or ant2 == @ant").times_amp.values))
        gains_len[ant]["phase"] = len(set(df_fitted.query("ant1 == @ant or ant2 == @ant").times_phase.values))
        gains_post[ant] = dict()
        gains_post[ant]["amp"] = dict()
        gains_post[ant]["phase"] = dict()
        gains_post[ant]["mean_phase"] = dict()
        gains_post[ant]["mean_phase"] = samples[:, j+gains_len[ant]["amp"]]
        for i, t in enumerate(df_fitted.query("ant1 == @ant or ant2 == @ant").times_amp.unique()):
            gains_post[ant]["amp"][t] = samples[:, j+i]
            # +1 mean skip ``mean_phase``
            if add_mean_phase:
                gains_post[ant]["phase"][t] = samples[:, j+gains_len[ant]["amp"]+i+1] + gains_post[ant]["mean_phase"]
            else:
                gains_post[ant]["phase"][t] = samples[:, j+gains_len[ant]["amp"]+i+1]

        j += gains_len[ant]["amp"] + gains_len[ant]["phase"] + 1

    fig, axes = plt.subplots(len(gains_len), 2, sharex=True, figsize=(8, 20))
    for i, ant in enumerate(gains_post):
        print("Plotting antenna ", ant)
        for t in gains_post[ant]["amp"].keys():
            # Arry with posterior for given t
            amp = gains_post[ant]["amp"][t]
            axes[i, 0].plot([t]*len(amp), amp, '.', color="#1f77b4")
        for t in gains_post[ant]["phase"].keys():
            phase = gains_post[ant]["phase"][t]
            axes[i, 1].plot([t]*len(phase), phase, '.', color="#1f77b4")
        axes[i, 1].yaxis.set_ticks_position("right")
    axes[0, 0].set_title("Amplitudes")
    axes[0, 1].set_title("Phases")
    axes[i, 0].set_xlabel("time, s")
    axes[i, 1].set_xlabel("time, s")
    if plotfn:
        fig.savefig(plotfn, bbox_inches="tight", dpi=100)
    fig.show()
    return fig


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
    # uvfits_fname = "/home/ilya/github/time_machine/bsc/reals/uvfs/J2038+5119_S_2005_07_20_yyk_uve.fits"
    # uvfits_fname = "/home/ilya/github/time_machine/bsc/reals/uvfs/BLLAC_RA_times.uvf"
    # uvfits_fname = "/home/ilya/github/time_machine/bsc/reals/uvfs/1502/1502_30s.uvf"
    uvfits_fname = "/home/ilya/github/time_machine/bsc/reals/J0005/J0005_30s.uvf"

    # for STOKES in range(2):
    #     for IF in range(2):
    STOKES = 0
    IF = 0
    # data_only_fname = "/home/ilya/github/time_machine/bsc/reals/RA/tests/BLLAC_STOKES_{}_IF_{}_amp120_phase60.txt".format(STOKES, IF)
    # data_only_fname = "/home/ilya/github/time_machine/bsc/reals/1502/1502_STOKES_{}_IF_{}_amp60_phase30.txt".format(STOKES, IF)
    data_only_fname = "/home/ilya/github/time_machine/bsc/reals/J0005/J0005_amp30_phase30_aver30.txt"
    # df = create_data_file(uvfits_fname, data_only_fname, STOKES=STOKES, IF=IF, step_amp=60, step_phase=30,
    #                       use_scans_for_amplitudes=False, calculate_noise=False)
                          # antennas_to_skip=(3, 8, 12, 13, 14, 16, 17))
    df = create_data_file_many_IFs(uvfits_fname, data_only_fname, step_amp=30, step_phase=30, calculate_noise=True)
    import sys
    sys.exit(0)
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
    re, im = gaussian_circ_ft(flux=2.0, dx=0.0, dy=0.0, bmaj=0.03, uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im
    re, im = gaussian_circ_ft(flux=1.0, dx=0.2, dy=0.0, bmaj=0.06, uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im
    re, im = gaussian_circ_ft(flux=0.5, dx=0.5, dy=0.25, bmaj=0.2, uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im

    # Plot
    fig = radplot(df, label="Sky Model")

    # Add gains
    df_updated, gains_dict = inject_gains(df,
                                          scale_gpamp=np.exp(7),
                                          scale_gpphase=np.exp(5),
                                          amp_gpamp=np.exp(-3),
                                          amp_gpphase=np.exp(-2))
    # Add noise
    # gains_noise_fname = "/home/ilya/github/time_machine/bsc/test.txt"
    gains_noise_fname = "/home/ilya/github/time_machine/bsc/test_int30_amp30_phase30.txt"
    df_updated = add_noise(df_updated, use_global_median_noise=True, use_per_baseline_median_noise=False,
                           outfname=gains_noise_fname)

    fig = radplot(df_updated, color="#ff7f0e", fig=fig, label="With gains")
    fig.savefig("J2001_check_models.png")
    fig = plot_model_gains(gains_dict)
    fig.savefig("J2001_check_gains.png")

    import json
    # Save gains
    with open("/home/ilya/github/time_machine/bsc/J2001_check_gains_30s.json", "w") as fo:
        json.dump(str(gains_dict), fo)
    # # Load gains
    # with open("/home/ilya/github/bsc/gains_0716.json", "r") as fo:
    #     loaded_gains_dict = json.load(fo)
    # import ast
    # loaded_gains_dict = ast.literal_eval(loaded_gains_dict)
