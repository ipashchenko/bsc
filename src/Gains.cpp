#include <iostream>
#include "Gains.h"
#include "Gain.h"
#include "MyExceptions.h"



Gains::Gains(Data data, int refant_) {
    // ant_i of the reference antenna
    refant = refant_;

    // Map between ant_i and its position in antennas vector.
    auto antennas_map = data.get_antennas_map();
    if(!antennas_map.count(refant)) {
        throw BadReferenceAntenna();
    }

    for (int i=0; i<data.n_antennas(); i++) {
        if (i != antennas_map[refant]) {
            antennas_changing_gain.push_back(i);
        }
    }

    // Map between antenna position in antennas vector and ant_i.
    auto antennas_map_inv = data.get_antennas_map_inv();
    for (int i=0; i<data.n_antennas(); i++) {
        gains.emplace_back(new Gain(antennas_map_inv[i], data.get_times_amp()[antennas_map_inv[i]], data.get_times_phase()[antennas_map_inv[i]]));
    }
}


Gains::Gains(Gains &other) {
    gains = std::vector<Gain*>();
    antennas_changing_gain = other.antennas_changing_gain;
    refant = other.refant;
    for (auto gain : other.gains) {
        gains.emplace_back(new Gain(*gain));
    }
}


Gains& Gains::operator=(const Gains& other) {
    if (&other == this)
        return *this;

    for (auto gain : gains) {
        delete gain;
    }

    gains = std::vector<Gain*>();
    for (auto gain : other.gains) {
        gains.emplace_back(new Gain(*gain));
    }
    antennas_changing_gain = other.antennas_changing_gain;
    refant = other.refant;
    return *this;
}


Gains::~Gains() {
    for (auto gain : gains) {
        delete gain;
    }
}


int Gains::size() const {
    return gains.size();
}


Gain* Gains::operator[](int i) {
    return gains[i];
}


void Gains::from_prior_amp_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_amp_mean(rng);
    }
}


void Gains::from_prior_phase_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_phase_mean(rng);
    }
}


void Gains::from_prior_v_amp(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_v_amp(rng);
    }
}


void Gains::from_prior_v_phase(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_v_phase(rng);
    }
}


void Gains::from_prior_hp_amp(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_hp_amp(rng);
    }
}


void Gains::from_prior_hp_phase(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_hp_phase(rng);
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
    int which = rng.rand_int(antennas_changing_gain.size());
    logH += gains[antennas_changing_gain[which]]->perturb(rng);
    return logH;
}


std::string Gains::description() const {
    std::string descr;
    for (auto gain : gains) {
        descr += gain->description();
        descr += " ";
    }
    descr.pop_back();
    return descr;
}


void Gains::print(std::ostream &out) const {
    for (auto gain : gains) {
        gain->print(out);
    }
}


