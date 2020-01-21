#include <iostream>
#include "DNestModel.h"
#include "RNG.h"


DNestModel::DNestModel() : logjitter(0.0), components(4, 30, false, MyConditionalPrior(30), DNest4::PriorType::log_uniform) {
    int refant_ant_i = 1;
    gains = new Gains(Data::get_instance(), refant_ant_i);
}


DNestModel::~DNestModel() {
    delete gains;
}


DNestModel::DNestModel(const DNestModel& other) : components(4, 30, false, MyConditionalPrior(30), DNest4::PriorType::log_uniform) {
    components = other.components;
    gains = new Gains(*other.gains);
    logjitter = other.logjitter;
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
    mu_real_full = other.mu_real_full;
    mu_imag_full = other.mu_imag_full;
}


DNestModel& DNestModel::operator=(const DNestModel& other) {
    if (this != &other) {
        components = other.components;
        *(gains) = *(other.gains);
        logjitter = other.logjitter;
        mu_real = other.mu_real;
        mu_imag = other.mu_imag;
        mu_real_full = other.mu_real_full;
        mu_imag_full = other.mu_imag_full;
    }
    return *this;
}


void DNestModel::from_prior(DNest4::RNG &rng) {
    // I do this because in ``calculate_sky_mu`` ``mu_real`` and ``mu_imag`` are multiplied and added.
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;
    mu_imag_full = zero;
    mu_real_full = zero;
    logjitter = -4.0 + 1.0*rng.randn();
    components.from_prior(rng);
    recenter();
    gains->from_prior_hp_amp(rng);
    gains->from_prior_hp_phase(rng);
    gains->from_prior_amp_mean(rng);
    gains->from_prior_phase_mean(rng);
    gains->from_prior_v_amp(rng);
    gains->from_prior_v_phase(rng);
    // Calculate C, L matrixes for rotation of v
    gains->calculate_C_amp();
    gains->calculate_C_phase();
    gains->calculate_L_amp();
    gains->calculate_L_phase();
    // Calculate latent v using calculated L and generated from prior v
    gains->calculate_amplitudes();
    gains->calculate_phases();
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

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else {
            logH = 0.0;
        }
        // No need to re-calculate sky_model or full model. Just calculate loglike.
        return logH;
    }

    // Perturb SkyModel
    else if(0.05 < u && u <= 0.3) {
        logH += components.perturb(rng);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        } else {
            logH = 0.0;
        }
        // This shouldn't be called in case of pre-rejection
        recenter();
        calculate_sky_mu(false);
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
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    std::valarray<double> theta;
    double c;
    std::valarray<double> b;
    std::valarray<double> ft;

    // Get the components
    const std::vector<std::vector<double>>& comps = (update)?(components.get_added()):
                                                    (components.get_components());

    if(!update)
    {
        // Zero the sky model prediction
        mu_real *= 0.0;
        mu_imag *= 0.0;
    }

    for(auto comp: comps)
    {
        // Phase shift due to not being in a phase center
        theta = 2*M_PI*mas_to_rad*(u*comp[0]+v*comp[1]);
        // Calculate FT of a Gaussian in a phase center
        c = pow(M_PI*exp(comp[3])*mas_to_rad, 2)/(4.*log(2.));
        ft = exp(comp[2]-c*(u*u + v*v));
        // Prediction of visibilities
        mu_real += ft*cos(theta);
        mu_imag += ft*sin(theta);
    }
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

    // Container for amplitudes and phases of gains for corresponding visibilities (TB) - means are inserted inside
    // individual Gain.calculate_amplitudes
    std::valarray<double> amp_ant_i (0.0, ant_i.size());
    std::valarray<double> amp_ant_j (0.0, ant_i.size());
    std::valarray<double> phase_ant_i (0.0, ant_i.size());
    std::valarray<double> phase_ant_j (0.0, ant_i.size());

    for (int k=0; k<ant_i.size(); k++) {
        auto ant_ik = antennas_map[ant_i[k]];
        auto ant_jk = antennas_map[ant_j[k]];
        auto idx_ik = idx_amp_ant_i[k];
        auto idx_jk = idx_amp_ant_j[k];
        auto ifk = IF[k];
        amp_ant_i[k] += gains->operator[](ant_ik)[ifk]->get_amplitudes()[idx_ik];
        amp_ant_j[k] += gains->operator[](ant_jk)[ifk]->get_amplitudes()[idx_jk];
        phase_ant_i[k] += gains->operator[](ant_ik)[ifk]->get_phases()[idx_ik];
        phase_ant_j[k] += gains->operator[](ant_jk)[ifk]->get_phases()[idx_jk];
    }

    // SkyModel modified by gains
    // This implies g = amp*exp(+1j*phi) and V = amp*exp(+1j*phi) - note "+" sign
    auto diff = phase_ant_i-phase_ant_j;
    auto cos_diff = cos(diff);
    auto sin_diff = sin(diff);
    auto amp_amp = amp_ant_i*amp_ant_j;
    mu_real_full = amp_amp*(cos_diff*mu_real - sin_diff*mu_imag);
    mu_imag_full = amp_amp*(cos_diff*mu_imag + sin_diff*mu_real);
}


double DNestModel::log_likelihood() const {
    // Grab the visibilities from the dataset
    const std::valarray<double>& vis_real = Data::get_instance().get_vis_real();
    const std::valarray<double>& vis_imag = Data::get_instance().get_vis_imag();
    const std::valarray<double>& sigma = Data::get_instance().get_sigma();

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
    components.print(out); out << '\t';
    gains->print(out);

}


std::string DNestModel::description() const
{
    std::string descr;

    // Anything printed by DNestModel::print (except the last line)
    descr += "logjitter ";

    // The rest is all what happens when you call .print on an RJObject
    descr += " dim_components max_num_components ";

    // Then the hyperparameters (i.e. whatever MyConditionalPrior::print prints)
    descr += " typical_flux dev_log_flux typical_radius dev_log_radius ";

    // Then the actual number of components
    descr += " num_components ";

    // Then it's all the components, padded with zeros
    // max_num_components is known in this model, so that's how far the
    // zero padding goes.
    for(int i=0; i<30; ++i)
        descr += " x[" + std::to_string(i) + "] ";
    for(int i=0; i<30; ++i)
        descr += " y[" + std::to_string(i) + "] ";
    for(int i=0; i<30; ++i)
        descr += " logflux[" + std::to_string(i) + "] ";
    for(int i=0; i<30; ++i)
        descr += " logbmaj[" + std::to_string(i) + "] ";

    descr += gains->description();

    return descr;
}



std::pair<double, double> DNestModel::center_mass() const {
    double x_b = 0;
    double y_b = 0;
    double flux_b = -100;
    for (auto comp : components.get_components()) {
        if(comp[2] > flux_b) {
            x_b = comp[0];
            y_b = comp[1];
            flux_b = comp[2];
        }
    }
    return std::make_pair<double, double>(reinterpret_cast<double &&>(x_b), reinterpret_cast<double &&>(y_b));
}


void DNestModel::shift_xy(std::pair<double, double> xy) {
    const std::vector<std::vector<double>>& comps = components.get_components();
    std::vector<std::vector<double>> new_comps;
    for (auto comp: comps) {
        comp[0] -= xy.first;
        comp[1] -= xy.second;
        new_comps.push_back(comp);
    }
    components.set_components(new_comps);
}


void DNestModel::recenter() {
    auto xy = center_mass();
    shift_xy(xy);
}