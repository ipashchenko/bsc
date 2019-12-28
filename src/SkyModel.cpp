#include <SkyModel.h>
#include <iostream>
#include "Data.h"

SkyModel::SkyModel() = default;


SkyModel::SkyModel(const SkyModel &other) {
    //std::cout << "In SkyModel copy ctor" << std::endl;
    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    components_sizes_ = other.components_sizes_;
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
    sin_theta = other.sin_theta;
    cos_theta = other.cos_theta;
}


SkyModel& SkyModel::operator=(const SkyModel& other) {
    //std::cout << "In SkyModel =" << std::endl;
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
    sin_theta = other.sin_theta;
    cos_theta = other.cos_theta;
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


void SkyModel::ft(const std::valarray<double>& u, const std::valarray<double>& v) {
    std::valarray<double> real (0.0, u.size());
    std::valarray<double> imag (0.0, u.size());
    // Traverse components and sum differences of new and old predictions for updated components.
    // TODO: Just choose updated component without cycle
    for (auto comp : components_) {
        if(comp->is_updated) {
            comp->ft(u, v);
            real = comp->get_mu_real() - comp->get_mu_real_old();
            imag = comp->get_mu_imag() - comp->get_mu_imag_old();
            if (comp->is_updated) {
                comp->is_updated = false;
                comp->update_old();
            }
            break;
        }
    }
    mu_real += real;
    mu_imag += imag;

}


void SkyModel::ft_from_all(const std::valarray<double>& u, const std::valarray<double>& v) {
    std::valarray<double> real (0.0, u.size());
    std::valarray<double> imag (0.0, u.size());
    for (auto comp : components_) {
        comp->ft(u, v);
        real = real + comp->get_mu_real();
        imag = imag + comp->get_mu_imag();
        comp->is_updated = false;
        comp->update_old();
    }
    mu_real = real;
    mu_imag = imag;
}


void SkyModel::set_param_vec(std::valarray<double> param) {
    size_t size_used = 0;

    for (auto comp : components_) {
        comp->set_param_vec(param[std::slice(size_used, comp->size(), 1)]);
        size_used += comp->size();
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
    sin_theta = zero;
    cos_theta = zero;
    for (auto comp: components_) {
        comp->from_prior(rng);
        comp->is_updated = true;
    }
    // Move center of mass to the phase center
    auto xy = center_brightest();
    shift_xy(xy);
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(get_n_components());
    components_[which]->is_updated = true;
    auto result = components_[which]->perturb(rng);
    update_brightest();
    return result;
}


int SkyModel::get_n_components() {
    return components_.size();
}


std::vector<int> SkyModel::get_components_sizes() {
    return components_sizes_;
}


void SkyModel::set_x(double x) {
    components_[0]->set_x(x);
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


void SkyModel::update_brightest() {
    double flux_b = 0;
    // Find brightest flux
    for (auto comp: components_) {
        if (comp->get_flux() > flux_b) {
            flux_b = comp->get_flux();
        }
    }
    // Updated labels
    for (auto comp: components_) {
        comp->is_brightest = comp->get_flux()==flux_b;
    }
}

std::pair<double, double> SkyModel::center_brightest() const {
    double x_b = 0;
    double y_b = 0;
    double flux_b = 0;
    for (auto comp: components_) {
        if(comp->get_flux() > flux_b) {
            x_b = comp->get_x();
            y_b = comp->get_y();
            flux_b = comp->get_flux();
        }
    }
    return std::make_pair<double, double>(reinterpret_cast<double &&>(x_b), reinterpret_cast<double &&>(y_b));
}


void SkyModel::shift_xy(std::pair<double, double> xy) {
    for (auto comp: components_) {
        comp->shift_xy(xy);
        // Phase shifting old component prediction because position of component relative to the phase center changed
        comp->phase_shift_old(cos_theta, sin_theta);
    }
}


void SkyModel::recenter_brightest() {
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    auto xy = center_brightest();
    auto theta = 2*M_PI*mas_to_rad*(-u*xy.first-v*xy.second);
    cos_theta = cos(theta);
    sin_theta = sin(theta);

    shift_xy(xy);
}


