#include <SkyModel.h>
#include <iostream>
#include <algorithm>
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
    for (auto comp : components_) {
        if(comp->is_updated) {
            comp->ft(u, v);
            real = comp->get_mu_real() - comp->get_mu_real_old();
            imag = comp->get_mu_imag() - comp->get_mu_imag_old();
            comp->is_updated = false;
            comp->update_old();
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
        if (comp->is_updated) {
            comp->is_updated = false;
            comp->update_old();
        }
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
    //std::cout << "Generating from prior SkyModel" << std::endl;
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;
    for (auto comp: components_) {
        comp->from_prior(rng);
        comp->is_updated = true;
    }
    // Move center of mass to the phase center
    auto xy = center_mass();
    shift_xy(xy);
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(get_n_components());
    components_[which]->is_updated = true;
    return components_[which]->perturb(rng);
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


std::pair<double, double> SkyModel::center_mass() const {
    double x = 0;
    double y = 0;
    double flux = 0;
    double sum_flux = 0;
    for (auto comp: components_) {
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
    auto xy = center_mass();
    shift_xy(xy);
}


bool SkyModel::are_sorted() {
    std::vector<double> distances;
    for (auto comp: components_) {
        auto x = comp->get_x();
        auto y = comp->get_y();
        distances.push_back(x*x + y*y);
    }
    return std::is_sorted(distances.begin(), distances.end());
}


