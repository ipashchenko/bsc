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

    gains = std::vector<std::vector<Gain*>>(data.n_antennas());
    // This map from antenna numbers to antenna indexes.
    auto antennas_map_inv = data.get_antennas_map_inv();
    for (int i=0; i<data.n_antennas(); i++) {
        for (int j=0; j<data.n_IF(); j++) {
            gains[i].emplace_back(new Gain(antennas_map_inv[i], j+1, data.get_times_amp()[antennas_map_inv[i]], data.get_times_phase()[antennas_map_inv[i]]));}
    }
}


Gains::Gains(Gains &other) {
    antennas_changing_gain = other.antennas_changing_gain;
    refant = other.refant;
    gains = std::vector<std::vector<Gain*>>(other.gains.size());
    for (int i=0; i < other.gains.size(); i++) {
        for (auto gain: other.gains[i]) {
            gains[i].emplace_back(new Gain(*gain));
        }
    }
}


Gains& Gains::operator=(const Gains& other) {
    if (&other == this)
        return *this;

    for (const auto & gain_vector : gains) {
        for (auto gain: gain_vector) {
            delete gain;
        }
    }

    for (int i=0; i < other.gains.size(); i++) {
        gains[i].resize(other.gains[i].size());
        for (int j=0; j < other.gains[i].size(); j++) {
            gains[i][j] = new Gain(*other.gains[i][j]);
        }
    }
    antennas_changing_gain = other.antennas_changing_gain;
    refant = other.refant;
    return *this;
}


Gains::~Gains() {
    for (const auto & gain_vector : gains) {
        for (auto gain: gain_vector) {
            delete gain;
        }
    }
}


int Gains::size() const {
    return gains.size();
}


std::vector<Gain*>& Gains::operator[](int i) {
    return gains[i];
}


void Gains::from_prior_amp_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        for (auto gain: gains[i]) {
            gain->from_prior_amp_mean(rng);
        }
    }
}


void Gains::from_prior_v_amp(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        for (auto gain: gains[i]) {
            gain->from_prior_v_amp(rng);
        }
    }
}


void Gains::from_prior_phase_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        for (auto gain: gains[i]) {
            gain->from_prior_phase_mean(rng);
        }
    }
}


void Gains::from_prior_v_phase(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        for (auto gain: gains[i]) {
            gain->from_prior_v_phase(rng);
        }
    }
}


void Gains::from_prior_hp_amp(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        for (auto gain: gains[i]) {
            gain->from_prior_hp_amp(rng);
        }
    }
}


void Gains::from_prior_hp_phase(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        for (auto gain: gains[i]) {
            gain->from_prior_v_amp(rng);
        }
    }
}


void Gains::calculate_C_amp() {
    for (const auto &gain_vector : gains) {
        for (auto gain: gain_vector) {
            gain->calculate_C_amp();
        }
    }
}


void Gains::calculate_C_phase() {
    for (const auto &gain_vector : gains) {
        for (auto gain: gain_vector) {
            gain->calculate_C_phase();
        }
    }
}


void Gains::calculate_L_amp() {
    for (const auto &gain_vector : gains) {
        for (auto gain: gain_vector) {
            gain->calculate_L_amp();
        }
    }
}


void Gains::calculate_L_phase() {
    for (const auto &gain_vector : gains) {
        for (auto gain: gain_vector) {
            gain->calculate_L_phase();
        }
    }
}


void Gains::calculate_amplitudes() {
    for (const auto &gain_vector : gains) {
        for (auto gain: gain_vector) {
            gain->calculate_amplitudes();
        }
    }
}


void Gains::calculate_phases() {
    for (const auto &gain_vector : gains) {
        for (auto gain: gain_vector) {
            gain->calculate_phases();
        }
    }
}


double Gains::perturb(DNest4::RNG& rng) {
    double logH = 0;
    int which_antenna = rng.rand_int(antennas_changing_gain.size());
    int which_IF = rng.rand_int(gains[which_antenna].size());
    logH += gains[antennas_changing_gain[which_antenna]][which_IF]->perturb(rng);

    return logH;
}


std::string Gains::description() const {
    std::string descr;
    for (const auto & gain_vector : gains) {
        for (auto gain: gain_vector) {
            descr += gain->description();
            descr += " ";
        }
    }
    descr.pop_back();
    return descr;
}


void Gains::print(std::ostream &out) const {
    for (const auto & gain_vector : gains) {
        for (auto gain: gain_vector) {
            gain->print(out);
        }
    }
}


