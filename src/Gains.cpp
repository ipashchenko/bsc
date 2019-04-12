#include <iostream>
#include "Gains.h"
#include "Gain.h"


Gains::Gains(Data data) {
    std::vector<int> antennas = data.get_instance().get_antennas();
    for (int i=0; i<data.n_antennas(); i++) {
        gains.emplace_back(new Gain(data.get_times_amp(), data.get_times_phase()));
    }
    std::cout << "Ctor of Gains finished" << std::endl;
}


void Gains::set_hp_amp(std::valarray<double> params) {
    size_t size_used = 0;
    for (auto gain : gains) {
        gain->set_hp_amp(params[std::slice(size_used, 2, 1)]);
        size_used += 2;
    }
}


void Gains::set_hp_phase(std::valarray<double> params) {
    size_t size_used = 0;
    for (auto gain : gains) {
        gain->set_hp_phase(params[std::slice(size_used, 2, 1)]);
        size_used += 2;
    }
}


void Gains::set_v_amp(std::valarray<double> params) {
    size_t size_used = 0;
    for (auto gain : gains) {
        gain->set_v_amp(params[std::slice(size_used, gain->size_amp(), 1)]);
        size_used += gain->size_amp();
    }
}


void Gains::set_v_phase(std::valarray<double> params) {
    size_t size_used = 0;
    for (auto gain : gains) {
        gain->set_v_phase(params[std::slice(size_used, gain->size_phase(), 1)]);
        size_used += gain->size_phase();
    }
}


int Gains::size() const {
    return gains.size();
}


Gain* Gains::operator[](int i) {
    return gains[i];


}


void Gains::from_prior_v_amp() {
    for (auto gain: gains) {
        gain->from_prior_v_amp();
    }
}


void Gains::from_prior_v_phase() {
    for (auto gain: gains) {
        gain->from_prior_v_phase();
    }
}


void Gains::from_prior_hp_amp(DNest4::RNG &rng) {
    std::cout << "Generating from prior Gains HP AMP" << std::endl;
    for (auto gain: gains) {
        gain->from_prior_hp_amp(rng);
    }
}


void Gains::from_prior_hp_phase(DNest4::RNG &rng) {
    for (auto gain: gains) {
        gain->from_prior_hp_phase(rng);
    }
}


void Gains::calculate_C_amp() {
    for (auto gain: gains) {
        gain->calculate_C_amp();
    }
}


void Gains::calculate_C_phase() {
    for (auto gain: gains) {
        gain->calculate_C_phase();
    }
}


void Gains::calculate_L_amp() {
    for (auto gain: gains) {
        gain->calculate_L_amp();
    }
}


void Gains::calculate_L_phase() {
    for (auto gain: gains) {
        gain->calculate_L_phase();
    }
}


void Gains::calculate_amplitudes() {
    for (auto gain: gains) {
        gain->calculate_amplitudes();
    }
}


void Gains::calculate_phases() {
    for (auto gain: gains) {
        gain->calculate_phases();
    }
}


void Gains::print_amplitudes(std::ostream &out) const {
    for (auto gain: gains) {
        gain->print_amplitudes(out);
    }
}


void Gains::print_phases(std::ostream &out) const {
    for (auto gain: gains) {
        gain->print_phases(out);
    }
}


void Gains::print_times(std::ostream &out) const {
    for (auto gain: gains) {
        gain->print_times(out);
    }
}


void Gains::print_hp(std::ostream &out) const {
    for (auto gain: gains) {
        gain->print_hp(out);
    }
}


void Gains::print_v(std::ostream &out) const {
    for (auto gain: gains) {
        gain->print_v(out);
    }
}


void Gains::print_C(std::ostream &out) const {
    for (auto gain: gains) {
        gain->print_C(out);
    }
}


void Gains::print_L(std::ostream &out) const {
    for (auto gain: gains) {
        gain->print_L(out);
    }
}


double Gains::perturb(DNest4::RNG& rng) {
    double logH = 0;
    int which = rng.rand_int(gains.size());
    logH += gains[which]->perturb(rng);
    std::cout << "Perturbing gain # " << which << " with logH =" << logH << std::endl;
    return logH;
}
