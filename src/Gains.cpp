#include <iostream>
#include "Gains.h"
#include "Gain.h"
#include "MyExceptions.h"
#include "Data.h"



Gains::Gains(Data& data, int refant_) {
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
    // FIXME
    gains.reserve(data.n_antennas());
    for (int i=0; i<data.n_antennas(); i++) {
        // FIXME
        gains.emplace_back(data.get_times_amp()[antennas_map_inv[i]], data.get_times_phase()[antennas_map_inv[i]]);
    }
}


int Gains::size() const {
    return gains.size();
}


Gain& Gains::operator[](int i) {
    return gains[i];
}


void Gains::from_prior_amp_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i].from_prior_amp_mean(rng);
    }
}


void Gains::from_prior_phase_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i].from_prior_phase_mean(rng);
    }
}


void Gains::from_prior_v_amp(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i].from_prior_v_amp(rng);
    }
}


void Gains::from_prior_v_phase(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i].from_prior_v_phase(rng);
    }
}


void Gains::from_prior_hp_amp(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i].from_prior_hp_amp(rng);
    }
}


void Gains::from_prior_hp_phase(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i].from_prior_hp_phase(rng);
    }
}


void Gains::calculate_C_amp() {
    for (auto& gain: gains) {
        gain.calculate_C_amp();
    }
}


void Gains::calculate_C_phase() {
    for (auto& gain: gains) {
        gain.calculate_C_phase();
    }
}


void Gains::calculate_L_amp() {
    for (auto& gain: gains) {
        gain.calculate_L_amp();
    }
}


void Gains::calculate_L_phase() {
    for (auto& gain: gains) {
        gain.calculate_L_phase();
    }
}


void Gains::calculate_amplitudes() {
    for (auto& gain: gains) {
        gain.calculate_amplitudes();
    }
}


void Gains::calculate_phases() {
    for (auto& gain: gains) {
        gain.calculate_phases();
    }
}


double Gains::perturb(DNest4::RNG& rng) {
    double logH = 0;
    int which = rng.rand_int(antennas_changing_gain.size());
    logH += gains[antennas_changing_gain[which]].perturb(rng);
    return logH;
}


std::string Gains::description() const {
    std::string descr;
    for (const auto& gain : gains) {
        descr += gain.description();
        descr += " ";
    }
    descr.pop_back();
    return descr;
}


void Gains::print(std::ostream &out) const {
    for (const auto& gain : gains) {
        gain.print(out);
    }
}
