#include <SkyModel.h>
#include <iostream>
#include "Data.h"


SkyModel::SkyModel(int n_comp) {
    for (int i=0; i<n_comp; i++) {
        components_.emplace_back(CGComponent());
    }
}

void SkyModel::ft(const std::valarray<double>& u, const std::valarray<double>& v) {
    // Zero prediction. However ``mu_real/imag`` must be already initialized. They are initialized to zeros in
    // ``SkyModel.from_prior`` and (in future) first ``ft_from_all`` should fill them with predictions.
    mu_real *= 0.0;
    mu_imag *= 0.0;
    for (auto& comp : components_) {
        comp.ft(u, v);
        mu_real += comp.get_mu_real();
        mu_imag += comp.get_mu_imag();
    }
}


void SkyModel::print(std::ostream &out) const
{
    for (const auto& comp: components_) {
        comp.print(out);
    }
}


void SkyModel::from_prior(DNest4::RNG &rng) {
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;

    for (auto& comp: components_) {
        comp.from_prior(rng);
    }
    // Move center of mass to the phase center
    auto xy = center_mass();
    shift_xy(xy);
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(components_.size());
    return components_[which].perturb(rng);
}


std::string SkyModel::description() const {
    std::string descr;
    for (const auto& comp: components_) {
        descr += comp.description();
        descr += " ";
    }
    descr.pop_back();
    return descr;
}


std::pair<double, double> SkyModel::center_mass() const {
    double x_b = 0;
    double y_b = 0;
    double flux_b = 0;
    for (const auto& comp: components_) {
        if(comp.get_flux() > flux_b) {
            x_b = comp.get_x();
            y_b = comp.get_y();
            flux_b = comp.get_flux();
        }
    }
    return std::make_pair<double, double>(reinterpret_cast<double &&>(x_b), reinterpret_cast<double &&>(y_b));
}


void SkyModel::shift_xy(std::pair<double, double> xy) {
    for (auto& comp: components_) {
        comp.shift_xy(xy);
    }
}


void SkyModel::recenter() {
    auto xy = center_mass();
    shift_xy(xy);
}
