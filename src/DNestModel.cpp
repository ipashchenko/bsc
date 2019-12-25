#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "RNG.h"
#include "Utils.h"


DNestModel::DNestModel() : logjitter(0.0) {

    sky_model = new SkyModel();
    int ncomp = 3;
    for (int i=0; i<ncomp; i++) {
        auto* comp = new CGComponent();
        sky_model->add_component(comp);
    }

    gains = new Gains(Data::get_instance());
    //std::cout << "gains # " << gains->size() << std::endl;
    //std::cout << "Ctor of DNestModel finished" << std::endl;
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
    //std::cout << "Generating from prior DNestModel" << std::endl;
    logjitter = -3.0 + 2.0*rng.randn();
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
    calculate_sky_mu(false);
    // Calculate full model (SkyModel + gains)
    calculate_mu();
    //std::cout << "End of DNestModel::from_prior - mu_real_full[0] = " << mu_real_full[0] << std::endl;
}


double DNestModel::perturb(DNest4::RNG &rng) {
    //std::cout << "In DNestModel::perturb" << std::endl;
    double logH = 0.;


    // Perturb jitter
    if(rng.rand() <= 0.1) {
        logH -= -0.5*pow((logjitter+3)/2.0, 2.0);
        logjitter += rng.randh();
        logH += -0.5*pow((logjitter+3)/2.0, 2.0);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;
        // No need to re-calculate sky_model or full model. Just calculate loglike.
    }

    // Perturb SkyModel
    if(rng.rand() <= 0.25) {
        logH += sky_model->perturb(rng);
        //std::cout << "Perturbed SkyModel with logH =" << logH << std::endl;

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            //std::cout << "Pre-rejected proposal SkyModel" << std::endl;
            return -1E300;
        }
        else
            logH = 0.0;
        // This shouldn't be called in case of pre-rejection
        sky_model->recenter();
        calculate_sky_mu(true);
    }
    // Perturb Gains
    else {
        //std::cout << "Perturbing gains" << std::endl;
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
    // It shouldn't be called in case of pre-rejection (if -1E300 is returned from sky_model or gains perturb
    calculate_mu();
    return logH;
}


void DNestModel::calculate_sky_mu(bool update) {
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    if(update and ft_calc_counter < 30) {
        sky_model->ft(u, v);
        ft_calc_counter += 1;
    } else {
        sky_model->ft_from_all(u, v);
        ft_calc_counter = 0;
    }
}


void DNestModel::calculate_mu() {

    // Mapping from antenna numbers (ant_i/j) to position in vector of antennas.
    std::unordered_map<int, int>& antennas_map = Data::get_instance().get_antennas_map();

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

    // Container for amplitudes and phases of gains for corresponding visibilities (TB) - means are inserted inside
    // individual Gain.calculate_amplitudes
    std::valarray<double> amp_ant_i (0.0, ant_i.size());
    std::valarray<double> amp_ant_j (0.0, ant_i.size());
    std::valarray<double> phase_ant_i (0.0, ant_i.size());
    std::valarray<double> phase_ant_j (0.0, ant_i.size());

    for (int k=0; k<ant_i.size(); k++) {
        amp_ant_i[k] += gains->operator[](antennas_map[ant_i[k]])->get_amplitudes()[idx_amp_ant_i[k]];
        amp_ant_j[k] += gains->operator[](antennas_map[ant_j[k]])->get_amplitudes()[idx_amp_ant_j[k]];
        phase_ant_i[k] += gains->operator[](antennas_map[ant_i[k]])->get_phases()[idx_phase_ant_i[k]];
        phase_ant_j[k] += gains->operator[](antennas_map[ant_j[k]])->get_phases()[idx_phase_ant_j[k]];
        //std::cout << "Amplitudes of gains in DNEstModel::calculate_mu : " << amp_ant_i[k] << ", " << amp_ant_j[k] << std::endl;
    }

    // SkyModel modified by gains
    // This implies g = amp*exp(+1j*phi) and V = amp*exp(+1j*phi) - note "+" sign
    mu_real_full = amp_ant_i*amp_ant_j*(cos(phase_ant_i-phase_ant_j)*sky_model_mu_real-
                                        sin(phase_ant_i-phase_ant_j)*sky_model_mu_imag);
    mu_imag_full = amp_ant_i*amp_ant_j*(cos(phase_ant_i-phase_ant_j)*sky_model_mu_imag+
                                        sin(phase_ant_i-phase_ant_j)*sky_model_mu_real);
}


double DNestModel::log_likelihood() const {
    //std::cout << "In DNestModel::log_likelihood" << std::endl;
    // Grab the visibilities from the dataset
    const std::valarray<double>& vis_real = Data::get_instance().get_vis_real();
    const std::valarray<double>& vis_imag = Data::get_instance().get_vis_imag();
    const std::valarray<double>& sigma = Data::get_instance().get_sigma();

    // Variance
    const std::valarray<double> var = sigma*sigma;
    //std::cout << "In DNestModel::log_likelihood - mu_real_full[0] = " << mu_real_full[0] << std::endl;
    // Complex Gaussian sampling distribution
    std::valarray<double> result = -log(2*M_PI*(var+exp(2.0*logjitter))) - 0.5*(pow(vis_real - mu_real_full, 2) +
        pow(vis_imag - mu_imag_full, 2))/(var+exp(2.0*logjitter))   ;
    double loglik = result.sum();
    //std::cout << "In DNestModel::log_likelihood - loglik = " << loglik << std::endl;
    return loglik;

}


void DNestModel::print(std::ostream &out) const {
    out << logjitter << '\t';
    sky_model->print(out);
    gains->print(out);
}


std::string DNestModel::description() const {
    std::string descr;
    descr += "logjitter ";
    descr += sky_model->description();
    descr += " ";
    descr += gains->description();

    return descr;
}


void DNestModel::set_x_skymodel(double x) {
    sky_model->set_x(x);
}

