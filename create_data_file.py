import numpy as np
import astropy.io.fits as pf
from astropy.time import Time
import pandas as pd


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


if __name__ == "__main__":
    uvfits_fname = "/home/ilya/github/DNest4/code/Examples/UV/J2001+2416_K_2006_06_11_yyk_vis.fits"
    df = create_data_file(uvfits_fname,
                          "/home/ilya/CLionProjects/ve2/uv_data_14.txt",
                          step_amp=120)
    print(df)