#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "RNG.h"
#include "Utils.h"


DNestModel::DNestModel() :
logjitter(0.0),
components(4, 10, false, MyConditionalPrior(-10, 10, -10, 10), DNest4::PriorType::log_uniform)
{

    //sky_model = new SkyModel();
    //int ncomp = 2;
    //for (int i=0; i<ncomp; i++) {
    //    auto* comp = new CGComponent();
    //    sky_model->add_component(comp);
    //}

    gains = new Gains(Data::get_instance());
    //std::cout << "gains # " << gains->size() << std::endl;
    //std::cout << "Ctor of DNestModel finished" << std::endl;
}


DNestModel::~DNestModel() {
    //delete sky_model;
    delete gains;
}



DNestModel::DNestModel(const DNestModel& other) : components(4, 10, false, MyConditionalPrior(-10, 10, -10, 10), DNest4::PriorType::log_uniform)
{
    //sky_model = new SkyModel(*other.sky_model);
    components = other.components;
    gains = new Gains(*other.gains);
    logjitter = other.logjitter;
    mu_real_full = other.mu_real_full;
    mu_imag_full = other.mu_imag_full;
}


DNestModel& DNestModel::operator=(const DNestModel& other) {
    if (this != &other) {
        //*(sky_model) = *(other.sky_model);
        components = other.components;
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
    //sky_model->from_prior(rng);
    components.from_prior(rng);
    recenter();
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
        // No need to re-calculate skymodel or full model. Just calculate loglike.
    }

    // Perturb SkyModel
    if(rng.rand() <= 0.25) {
        logH += components.perturb(rng);

        // Pre-reject
        if(rng.rand() >= exp(logH))
            return -1E300;
        else
            logH = 0.0;

        // This shouldn't be called in case of pre-rejection
        calculate_sky_mu(components.get_removed().size() == 0);
        recenter();
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

    std::valarray<double> theta;
    double c;
    std::valarray<double> b;
    std::valarray<double> ft;

    const std::vector<std::vector<double>>& comps = components.get_components();

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
    const std::vector<int>& idx_amp_ant_i = Data::get_instance().get_idx_amp_ant_i();
    const std::vector<int>& idx_amp_ant_j = Data::get_instance().get_idx_amp_ant_j();
    const std::vector<int>& idx_phase_ant_i = Data::get_instance().get_idx_phase_ant_i();
    const std::vector<int>& idx_phase_ant_j = Data::get_instance().get_idx_phase_ant_j();

    //// SkyModel predictions
    //std::valarray<double> sky_model_mu_real = sky_model->get_mu_real();
    //std::valarray<double> sky_model_mu_imag = sky_model->get_mu_imag();

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
    mu_real_full = amp_ant_i*amp_ant_j*(cos(phase_ant_i-phase_ant_j)*mu_real-
                                        sin(phase_ant_i-phase_ant_j)*mu_imag);
    mu_imag_full = amp_ant_i*amp_ant_j*(cos(phase_ant_i-phase_ant_j)*mu_imag+
                                        sin(phase_ant_i-phase_ant_j)*mu_real);
}


double DNestModel::log_likelihood() const {
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
    components.print(out); out<<'\t';
    gains->print(out);
}


std::string DNestModel::description() const {
    std::string descr;
    descr += "logjitter ";

    const std::vector<std::vector<double>>& comps = components.get_components();
    for (auto comp: comps) {
        descr += "x y logflux logbmaj ";
    }

    descr += gains->description();

    return descr;
}


std::pair<double, double> DNestModel::center_mass() const {
    double x = 0;
    double y = 0;
    double flux = 0;
    double sum_flux = 0;

    const std::vector<std::vector<double>>& comps = components.get_components();
    for(auto comp: comps)
    {
        flux = exp(comp[2]);
        x += comp[0]*flux;
        y += comp[1]*flux;
        sum_flux += flux;
    }
    return std::make_pair<double, double>(x/sum_flux, y/sum_flux);
}


void DNestModel::shift_xy(std::pair<double, double> xy) {
    const std::vector<std::vector<double>>& comps = components.get_components();
    for (auto comp: comps) {
        comp[0] -= xy.first;
        comp[1] -= xy.second;
    }
}


void DNestModel::recenter() {
    auto xy = center_mass();
    shift_xy(xy);
}