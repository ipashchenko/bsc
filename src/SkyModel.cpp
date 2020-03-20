#include <SkyModel.h>
#include <iostream>
#include "Data.h"

SkyModel::SkyModel() = default;


SkyModel::SkyModel(const SkyModel &other) {
    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    components_sizes_ = other.components_sizes_;
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
    idx_brightest = other.idx_brightest;
    idx_brightest_old = other.idx_brightest_old;
}


SkyModel& SkyModel::operator=(const SkyModel& other) {
    if (&other == this) {
        return *this;
    }

    for (auto comp : components_) {
        delete comp;
    }
    components_.clear();

    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    components_sizes_ = other.components_sizes_;
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
    idx_brightest = other.idx_brightest;
    idx_brightest_old = other.idx_brightest_old;
    return *this;
}



SkyModel::~SkyModel() {
    for (auto comp : components_) {
        delete comp;
    }
    components_.clear();
}


void SkyModel::add_component(Component *component) {
    components_.push_back(component);
    components_sizes_.push_back(component->size());
}


void SkyModel::phase_shift_mu(std::pair<double, double> shift) {
    //std::cout << "In SkyModel.phase_shift_mu ";
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    auto theta = 2*M_PI*mas_to_rad*(-u*shift.first-v*shift.second);
    auto cos_theta = cos(theta);
    auto sin_theta = sin(theta);
    mu_real = cos_theta*mu_real - sin_theta*mu_imag;
    mu_imag = cos_theta*mu_imag + sin_theta*mu_real;
}

void SkyModel::ft_update(const std::valarray<double>& u, const std::valarray<double>& v) {
    // All component are already shifted => phase shift old visibilities
    // Traverse components and sum differences of new and old predictions for updated components.
    //std::cout << "idx brightest = " << idx_brightest << ", old = " << idx_brightest_old << std::endl;
    for (auto comp : components_) {
        if(comp->is_updated) {
            // After possible shift (if needed)
            comp->ft(u, v);
            comp->update_old();
            // New brightest component OR position of the brightest component has changed
            if(idx_brightest =! idx_brightest_old || components_[idx_brightest]->is_position_updated) {
                //std::cout << "Last shift = " << last_shift.first << ", " << last_shift.second << std::endl;
                // Phase shift current sky model prediction
                phase_shift_mu(last_shift);
                // Phase shift old prediction of component;
                comp->phase_shift_mu_old(last_shift);
            }
            // TODO: Here I can add to ``mu_real``, ``mu_imag`` w/o declaring new arrays. However ``mu_real/imag`` must
            // be already initialized. They are initialized to zeros in ``SkyModel.from_prior`` and first
            // ``ft_from_all`` fills them with predictions.
            mu_real += (comp->get_mu_real() - comp->get_mu_real_old());
            mu_imag += (comp->get_mu_imag() - comp->get_mu_imag_old());
            comp->is_updated = false;
            comp->is_position_updated = false;
            //comp->update_old();
            break;
        }
    }
}


void SkyModel::ft(const std::valarray<double>& u, const std::valarray<double>& v) {
    // Zero prediction. However ``mu_real/imag`` must be already initialized. They are initialized to zeros in
    // ``SkyModel.from_prior`` and (in future) first ``ft_from_all`` should fill them with predictions.
    mu_real *= 0.0;
    mu_imag *= 0.0;
    for (auto comp : components_) {
        comp->ft(u, v);
        // TODO: Here I can add to ``mu_real``, ``mu_imag`` w/o declaring new arrays.
        mu_real += comp->get_mu_real();
        mu_imag += comp->get_mu_imag();
        comp->update_old();
    }
}


size_t SkyModel::size() const {
    size_t cursize = 0;
    if (components_.empty()) {
        return cursize;
    }
    else {
        for (auto comp : components_) {
            cursize += comp->size();
        }
    }
    return cursize;
}


void SkyModel::print(std::ostream &out) const
{
    for (auto comp: components_) {
        comp->print(out);
    }
}


void SkyModel::from_prior(DNest4::RNG &rng) {
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;

    for (auto comp: components_) {
        comp->from_prior(rng);
        comp->is_updated = true;
    }
    // Find brightest, move brightest component to the phase center if necessary
    recenter();
    //std::cout << "In SkyModel::from_prior brightest = " << idx_brightest << ", old = " << idx_brightest_old << std::endl;
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(get_n_components());
    components_[which]->is_updated = true;
    return components_[which]->perturb(rng);
}


int SkyModel::get_n_components() {
    return components_.size();
}



std::string SkyModel::description() const {
    std::string descr;
    for (auto comp: components_) {
        descr += comp->description();
        descr += " ";
    }
    descr.pop_back();
    return descr;
}


void SkyModel::find_brightest() {
    //std::cout << "Finding brightest component" << std::endl;
    if(idx_brightest != -1) {
        idx_brightest_old = idx_brightest;
    }
    double flux_b = 0;
    for(int i=0; i<components_.size(); i++) {
        // get_flux() returns exp of log_flux
        if(components_[i]->get_flux() > flux_b) {
            flux_b = components_[i]->get_flux();
            idx_brightest = i;
        }
    }
}

std::pair<double, double> SkyModel::center_mass() const {
    auto comp = components_[idx_brightest];
    return std::make_pair<double, double>(comp->get_x(), comp->get_y());
}


std::pair<double, double> SkyModel::center_mass2() const {
    double x = 0;
    double y = 0;
    double flux = 0;
    double sum_flux = 0;
    for (auto comp: components_) {
        // get_flux() returns exp of log_flux
        flux = comp->get_flux();
        x += comp->get_x()*flux;
        y += comp->get_y()*flux;
        sum_flux += flux;

    }
    return std::make_pair<double, double>(x/sum_flux, y/sum_flux);
}


void SkyModel::shift_xy(std::pair<double, double> xy) {
    for (auto comp: components_) {
        comp->shift_xy(xy);
    }
}


void SkyModel::recenter() {
    find_brightest();
    // Shift phase center only if the brightest component has changed OR its position has changd
    if (idx_brightest != idx_brightest_old || components_[idx_brightest]->is_position_updated) {
        auto xy = center_mass();
        shift_xy(xy);
        last_shift = xy;
    }
}


