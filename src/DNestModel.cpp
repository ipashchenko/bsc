#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "RNG.h"
#include "Utils.h"


DNestModel::DNestModel() {

    sky_model = new SkyModel();
    int ncomp = 3;
    for (int i=0; i<ncomp; i++) {
        auto* comp = new CGComponent();
        sky_model->add_component(comp);
    }

    gains = new Gains(Data::get_instance());
}


DNestModel::~DNestModel() {
    delete sky_model;
    delete gains;
}



DNestModel::DNestModel(const DNestModel& other) {
    sky_model = new SkyModel(*other.sky_model);
    gains = new Gains(*other.gains);
    logjitter = other.logjitter;
    mu_real_full = other.mu_real_full;
    mu_imag_full = other.mu_imag_full;
}


DNestModel& DNestModel::operator=(const DNestModel& other) {
    if (this != &other) {
        *(sky_model) = *(other.sky_model);
        *(gains) = *(other.gains);
        logjitter = other.logjitter;
        mu_real_full = other.mu_real_full;
        mu_imag_full = other.mu_imag_full;
    }
    return *this;
}


void DNestModel::from_prior(DNest4::RNG &rng) {
    for (int i=0;i<Data::get_instance().n_antennas();i++) {
        logjitter.push_back(-1.5 + 1.0*rng.randn());
    }
    sky_model->from_prior(rng);
    //sky_model->print(std::cout);
    gains->from_prior_hp_amp(rng);
    gains->from_prior_hp_phase(rng);
    //gains->print_hp(std::cout);
    gains->from_prior_v_amp();
    gains->from_prior_v_phase();
    //gains->print_v(std::cout);
    // Calculate C, L matrixes for rotation of v
    gains->calculate_C_amp();
    gains->calculate_C_phase();
    gains->calculate_L_amp();
    gains->calculate_L_phase();
    // Calculate latent v using calculated L and generated from prior v
    gains->calculate_amplitudes();
    gains->calculate_phases();
    //gains->print_amplitudes(std::cout);
    //gains->print_phases(std::cout);
    // Calculate SkyModel
    calculate_sky_mu();
    // Calculate full model (SkyModel + gains)
    calculate_mu();
}


double DNestModel::perturb(DNest4::RNG &rng) {
    double logH = 0.;

    double u = rng.rand();

    // Perturb jitter
    if(u <= 0.05) {
        int which_jitter = rng.rand_int(Data::get_instance().n_antennas());
        logH -= -0.5*pow((logjitter[which_jitter]+1.5)/1.0, 2.0);
        logjitter[which_jitter] += 1.0*rng.randh();
        logH += -0.5*pow((logjitter[which_jitter]+1.5)/1.0, 2.0);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;
        // No need to re-calculate sky_model or full model. Just calculate re-calculate loglike with new jitter.
        return logH;
    }

    // Perturb SkyModel
    else if(0.05 < u <= 0.25) {
        logH += sky_model->perturb(rng);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            //std::cout << "Pre-rejected proposal SkyModel" << std::endl;
            return -1E300;
        }
        else
            logH = 0.0;
        // This shouldn't be called in case of pre-rejection
        sky_model->recenter();
        calculate_sky_mu();
    }
    // Perturb Gains
    else {
        // C, L and v are re-calculated by individual Gains that is perturbed
        logH += gains->perturb(rng);
        // Gains pre-reject also in individual Gain instances
        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            //std::cout << "Pre-rejected proposal SkyModel" << std::endl;
            return -1E300;
        }
        else
            logH = 0.0;
    }
    // It shouldn't be called in case of pre-rejection (if -1E300 is returned from SkyModel or Gains perturb)
    calculate_mu();
    return logH;
}


void DNestModel::calculate_sky_mu() {
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();
    sky_model->ft(u, v);
}


void DNestModel::calculate_mu() {

    // Mapping from antenna numbers (ant_i/j) to position in vector of antennas.
    std::unordered_map<int, int>& antennas_map = Data::get_instance().get_antennas_map();

    // Grab antenna numbers and indexes in gain curves of individual visibility measurements
    const std::vector<int>& ant_i = Data::get_instance().get_ant_i();
    const std::vector<int>& ant_j = Data::get_instance().get_ant_j();
    const std::vector<int>& stokes = Data::get_instance().get_stokes();
    const std::vector<int>& IF = Data::get_instance().get_IF();
    const std::vector<int>& idx_amp_ant_i = Data::get_instance().get_idx_amp_ant_i();
    const std::vector<int>& idx_amp_ant_j = Data::get_instance().get_idx_amp_ant_j();
    const std::vector<int>& idx_phase_ant_i = Data::get_instance().get_idx_phase_ant_i();
    const std::vector<int>& idx_phase_ant_j = Data::get_instance().get_idx_phase_ant_j();

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
        amp_ant_i[k] += gains->operator[](antennas_map[ant_i[k]])[IF[k]]->get_amplitudes()[idx_amp_ant_i[k]];
        amp_ant_j[k] += gains->operator[](antennas_map[ant_j[k]])[IF[k]]->get_amplitudes()[idx_amp_ant_j[k]];
        phase_ant_i[k] += gains->operator[](antennas_map[ant_i[k]])[IF[k]]->get_phases()[idx_phase_ant_i[k]];
        phase_ant_j[k] += gains->operator[](antennas_map[ant_j[k]])[IF[k]]->get_phases()[idx_phase_ant_j[k]];
    }

    // SkyModel modified by gains
    // This implies g = amp*exp(+1j*phi) and V = amp*exp(+1j*phi) - note "+" sign
    mu_real_full = amp_ant_i*amp_ant_j*(cos(phase_ant_i-phase_ant_j)*sky_model_mu_real-
                                        sin(phase_ant_i-phase_ant_j)*sky_model_mu_imag);
    mu_imag_full = amp_ant_i*amp_ant_j*(cos(phase_ant_i-phase_ant_j)*sky_model_mu_imag+
                                        sin(phase_ant_i-phase_ant_j)*sky_model_mu_real);
}


double DNestModel::log_likelihood() const {
    // Grab the visibilities from the dataset
    const std::valarray<double>& vis_real = Data::get_instance().get_vis_real();
    const std::valarray<double>& vis_imag = Data::get_instance().get_vis_imag();
    const std::valarray<double>& sigma = Data::get_instance().get_sigma();

    const std::vector<int>& ant_i = Data::get_instance().get_ant_i();
    const std::vector<int>& ant_j = Data::get_instance().get_ant_j();
    std::valarray<double> logjitter_array (0.0, vis_real.size());
    for (int i=0; i<vis_real.size(); i++) {
        logjitter_array[i] = logjitter[ant_i[i]]+logjitter[ant_j[i]];
    }

    // Variance
    const std::valarray<double> var = sigma*sigma;
    // Complex Gaussian sampling distribution
    std::valarray<double> result = -log(2*M_PI*(var+exp(2.0*logjitter_array))) - 0.5*(pow(vis_real - mu_real_full, 2) +
        pow(vis_imag - mu_imag_full, 2))/(var+exp(2.0*logjitter_array))   ;
    double loglik = result.sum();
    return loglik;

}


void DNestModel::print(std::ostream &out) const {
    for (int i = 0; i < logjitter.size(); i++) {
        out << logjitter[i] << '\t';
    }
    sky_model->print(out);
    gains->print(out);
}


std::string DNestModel::description() const {
    std::string descr;
    for (int i = 0; i < logjitter.size(); i++) {
        descr += ("logjitter" + std::to_string(i) + " ");
    }
    descr += sky_model->description();
    descr += " ";
    descr += gains->description();

    return descr;
}
