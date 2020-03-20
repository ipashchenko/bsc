#include <iostream>
#include "DNestModel.h"
#include "RNG.h"


DNestModel::DNestModel() : logjitter(0.0), ft_calc_counter(0) {

    sky_model = new SkyModel();
    int ncomp = 7;
    for (int i=0; i<ncomp; i++) {
        auto* comp = new CGComponent();
        sky_model->add_component(comp);
    }

    int refant_ant_i = 1;
    gains = new Gains(Data::get_instance(), refant_ant_i);


    // Mapping from antenna numbers (ant_i/j) to position in vector of antennas.
    std::unordered_map<int, int>& antennas_map = Data::get_instance().get_antennas_map();

    const std::vector<int>& idx_amp_ant_i = Data::get_instance().get_idx_amp_ant_i();
    const std::vector<int>& idx_amp_ant_j = Data::get_instance().get_idx_amp_ant_j();
    const std::vector<int>& idx_phase_ant_i = Data::get_instance().get_idx_phase_ant_i();
    const std::vector<int>& idx_phase_ant_j = Data::get_instance().get_idx_phase_ant_j();

    const std::vector<int>& ant_i = Data::get_instance().get_ant_i();
    const std::vector<int>& ant_j = Data::get_instance().get_ant_j();

    for (int k=0; k<ant_i.size(); k++) {
        ant_ik.emplace_back(antennas_map[ant_i[k]]);
        ant_jk.emplace_back(antennas_map[ant_j[k]]);

        idx_ik_amp.emplace_back(idx_amp_ant_i[k]);
        idx_jk_amp.emplace_back(idx_amp_ant_j[k]);
        idx_ik_phase.emplace_back(idx_phase_ant_i[k]);
        idx_jk_phase.emplace_back(idx_phase_ant_j[k]);
    }
}


DNestModel::~DNestModel() {
    delete sky_model;
    delete gains;
}



DNestModel::DNestModel(const DNestModel& other) {
    sky_model = new SkyModel(*other.sky_model);
    gains = new Gains(*other.gains);
    logjitter = other.logjitter;
    ft_calc_counter = other.ft_calc_counter;
    mu_real_full = other.mu_real_full;
    mu_imag_full = other.mu_imag_full;
    ant_ik = other.ant_ik;
    ant_jk = other.ant_jk;
    idx_ik_amp = other.idx_ik_amp;
    idx_jk_amp = other.idx_jk_amp;
    idx_ik_phase = other.idx_ik_phase;
    idx_jk_phase = other.idx_jk_phase;
}


DNestModel& DNestModel::operator=(const DNestModel& other) {
    if (this != &other) {
        *(sky_model) = *(other.sky_model);
        *(gains) = *(other.gains);
        logjitter = other.logjitter;
        ft_calc_counter = other.ft_calc_counter;
        mu_real_full = other.mu_real_full;
        mu_imag_full = other.mu_imag_full;
        ant_ik = other.ant_ik;
        ant_jk = other.ant_jk;
        idx_ik_amp = other.idx_ik_amp;
        idx_jk_amp = other.idx_jk_amp;
        idx_ik_phase = other.idx_ik_phase;
        idx_jk_phase = other.idx_jk_phase;
    }
    return *this;
}


void DNestModel::from_prior(DNest4::RNG &rng) {
    //std::cout << "From prior DNestModel" << std::endl;
    logjitter = -4.0 + 1.0*rng.randn();
    //for (int i=0;i<Data::get_instance().n_antennas();i++) {
    //    logjitter.push_back(-2.0 + 0.5*rng.randn());
    //}
    sky_model->from_prior(rng);
    gains->from_prior_amp_mean(rng);
    gains->from_prior_phase_mean(rng);
    gains->from_prior_amp(rng);
    gains->from_prior_phase(rng);
    gains->sum();
    // Calculate SkyModel
    calculate_sky_mu(false);
    // Calculate full model (SkyModel + gains)
    calculate_mu();

}


double DNestModel::perturb(DNest4::RNG &rng) {
    double logH = 0.;

    double u = rng.rand();

    // Perturb jitter
    if(u <= 0.05) {
        logH -= -0.5*pow((logjitter+4)/1.0, 2.0);
        logjitter += 1.0*rng.randh();
        logH += -0.5*pow((logjitter+4)/1.0, 2.0);
        //int which_jitter = rng.rand_int(Data::get_instance().n_antennas());
        //logH -= -0.5*pow((logjitter[which_jitter]+2.0)/0.5, 2.0);
        //logjitter[which_jitter] += 0.5*rng.randh();
        //logH += -0.5*pow((logjitter[which_jitter]+2.0)/0.5, 2.0);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;
        // No need to re-calculate sky_model or full model. Just calculate loglike.
        return logH;
    }

    // Perturb SkyModel
    else if(0.05 < u && u <= 0.3) {
        logH += sky_model->perturb(rng);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        } else {
            logH = 0.0;
        }
        // This shouldn't be called in case of pre-rejection
        sky_model->recenter();
        calculate_sky_mu(true);
    }
    // Perturb Gains
    else {
        // C, L and v are re-calculated by individual Gains that is perturbed
        logH += gains->perturb(rng);
        // Gains pre-reject also in individual Gain instances
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        } else {
            logH = 0.0;
        }
    }
    // It shouldn't be called in case of pre-rejection (if -1E300 is returned from sky_model or gains perturb
    calculate_mu();
    return logH;
}


void DNestModel::calculate_sky_mu(bool update) {
    // Get (u, v)-points for calculating SkyModel predictions
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();
    // FT (calculate SkyModel prediction)
    if(update && ft_calc_counter < 30) {
        //std::cout << "Update" << std::endl;
        sky_model->ft_update(u, v);
        ft_calc_counter += 1;
    } else {
        //std::cout << "Full" << std::endl;
        sky_model->ft(u, v);
        ft_calc_counter = 0;
    }
    //std::cout << "COUNTER = " << ft_calc_counter << std::endl;
}


void DNestModel::calculate_mu() {

    // Grab antenna numbers and indexes in gain curves of individual visibility measurements
    const std::vector<int>& ant_i = Data::get_instance().get_ant_i();

    // SkyModel predictions
    std::valarray<double> sky_model_mu_real = sky_model->get_mu_real();
    std::valarray<double> sky_model_mu_imag = sky_model->get_mu_imag();

    // Container for amplitudes and phases of gains for corresponding visibilities (TB) - means are inserted inside
    // individual Gain.calculate_amplitudes
    std::valarray<double> amp_ant_i (0.0, ant_i.size());
    std::valarray<double> amp_ant_j (0.0, ant_i.size());
    std::valarray<double> phase_ant_i (0.0, ant_i.size());
    std::valarray<double> phase_ant_j (0.0, ant_i.size());


    for (int k=0; k<ant_i.size(); k++) {
        amp_ant_i[k] += gains->operator[](ant_ik[k])->get_amplitudes()[idx_ik_amp[k]];
        amp_ant_j[k] += gains->operator[](ant_jk[k])->get_amplitudes()[idx_jk_amp[k]];
        phase_ant_i[k] += gains->operator[](ant_ik[k])->get_phases()[idx_ik_phase[k]];
        phase_ant_j[k] += gains->operator[](ant_jk[k])->get_phases()[idx_jk_phase[k]];

    }

    // SkyModel modified by gains
    // This implies g = amp*exp(+1j*phi) and V = amp*exp(+1j*phi) - note "+" sign
    auto diff = phase_ant_i-phase_ant_j;
    auto cos_diff = cos(diff);
    auto sin_diff = sin(diff);
    auto amp_amp = amp_ant_i*amp_ant_j;
    mu_real_full = amp_amp*(cos_diff*sky_model_mu_real - sin_diff*sky_model_mu_imag);
    mu_imag_full = amp_amp*(cos_diff*sky_model_mu_imag + sin_diff*sky_model_mu_real);
}


double DNestModel::log_likelihood() const {
    // Grab the visibilities from the dataset
    const std::valarray<double>& vis_real = Data::get_instance().get_vis_real();
    const std::valarray<double>& vis_imag = Data::get_instance().get_vis_imag();
    const std::valarray<double>& sigma = Data::get_instance().get_sigma();

    //const std::vector<int>& ant_i = Data::get_instance().get_ant_i();
    //const std::vector<int>& ant_j = Data::get_instance().get_ant_j();
    //std::valarray<double> logjitter_array (0.0, vis_real.size());
    //for (int i=0; i<vis_real.size(); i++) {
    //    logjitter_array[i] = logjitter[ant_i[i]]+logjitter[ant_j[i]];
    //}

    // Variance
    const std::valarray<double> var = sigma*sigma;
    // Complex Gaussian sampling distribution
    std::valarray<double> result = -log(2*M_PI*(var+exp(2.0*logjitter))) - 0.5*(pow(vis_real - mu_real_full, 2) +
        pow(vis_imag - mu_imag_full, 2))/(var+exp(2.0*logjitter))   ;
    double loglik = result.sum();
    return loglik;

}


void DNestModel::print(std::ostream &out) const {
    out << logjitter << '\t';
    //for (double i : logjitter) {
    //    out << i << '\t';
    //}
    sky_model->print(out);
    gains->print(out);
}


std::string DNestModel::description() const {
    std::string descr;
    descr += "logjitter ";
    //auto antennas_map_inv = Data::get_instance().get_antennas_map_inv();
    //// logjitter array must have the same length as number of antennas
    //for (int i = 0; i < Data::get_instance().n_antennas(); i++) {
    //    descr += ("logjitter_" + std::to_string(antennas_map_inv[i]) + " ");
    //}
    descr += sky_model->description();
    descr += " ";
    descr += gains->description();

    return descr;
}

