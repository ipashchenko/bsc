#include "DNestModel.h"
#include "RNG.h"


DNestModel::DNestModel() {

    sky_model = new SkyModel();
    int ncomp = 1;
    for (int i=0; i<ncomp; i++) {
        auto* comp = new CGComponent();
        sky_model->add_component(comp);
    }

    gains = new Gains(Data::get_instance());

}


void DNestModel::from_prior(DNest4::RNG &rng) {
    sky_model->from_prior(rng);
    gains->from_prior_hp_amp(rng);
    gains->from_prior_hp_phase(rng);
    gains->from_prior_v_amp();
    gains->from_prior_v_phase();
    // Calculate C, L matrixes for rotation of v
    gains->calculate_C_amp();
    gains->calculate_C_phase();
    gains->calculate_L_amp();
    gains->calculate_L_phase();
    // Calculate latent v using calculated L and generated from prior v
    gains->calculate_amplitudes();
    gains->calculate_phases();
    // Calculate SkyModel
    calculate_sky_mu();
    // Calculate full model (SkyModel + gains)
    calculate_mu();
}


double DNestModel::perturb(DNest4::RNG &rng) {
    double logH = 0.;
    // Perturb SkyModel
    if(rng.rand() <= 0.25) {
        logH += sky_model->perturb(rng);
    }
    // Perturb Gains
    else {
        logH += gains->perturb(rng);
    }
    return logH;
}


void DNestModel::calculate_sky_mu() {
    // Get (u, v)-points for calculating SkyModel predictions
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();
    // FT (calculate SkyModel prediction)
    sky_model->ft(u, v);
}


void DNestModel::calculate_mu() {

    // Grab antenna numbers and indexes in gain curves of individual visibility measurements
    const std::vector<int>& ant_i = Data::get_instance().get_ant_i();
    const std::vector<int>& ant_j = Data::get_instance().get_ant_j();
    const std::vector<int>& idx_amp_ant_i = Data::get_instance().get_idx_amp_ant_i();
    const std::vector<int>& idx_amp_ant_j = Data::get_instance().get_idx_amp_ant_j();
    const std::vector<int>& idx_phase_ant_i = Data::get_instance().get_idx_phase_ant_i();
    const std::vector<int>& idx_phase_ant_j = Data::get_instance().get_idx_phase_ant_j();

    // SkyModel predictions
    std::valarray<double> sky_model_mu_real = sky_model->get_mu_real();
    std::valarray<double> sky_model_mu_imag = sky_model->get_mu_imag();

    // Container for amplitudes and phases of gains for corresponding visibilities (TB) - means are already set up
    std::valarray<double> amp_ant_i (1.0, ant_i.size());
    std::valarray<double> amp_ant_j (1.0, ant_i.size());
    std::valarray<double> phase_ant_i (0.0, ant_i.size());
    std::valarray<double> phase_ant_j (0.0, ant_i.size());

    for (int k=0; k<ant_i.size(); k++) {
        amp_ant_i[k] += gains->operator[](ant_i[k]-1)->get_amplitudes()[idx_amp_ant_i[k]];
        amp_ant_j[k] += gains->operator[](ant_j[k]-1)->get_amplitudes()[idx_amp_ant_j[k]];
        phase_ant_i[k] += gains->operator[](ant_i[k]-1)->get_phases()[idx_amp_ant_i[k]];
        phase_ant_j[k] += gains->operator[](ant_j[k]-1)->get_phases()[idx_amp_ant_j[k]];
    }

    // SkyModel modified by gains
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

    // Variance
    const std::valarray<double> var = sigma*sigma;

    // Complex Gaussian sampling distribution
    std::valarray<double> result = -log(2*PI*var) - 0.5*(pow(vis_real - mu_real_full, 2) +
        pow(vis_imag - mu_imag_full, 2))/var;
    return result.sum();

}


void DNestModel::print(std::ostream &out) const {

}


std::string DNestModel::description() const {
    return std::__cxx11::string();
}
