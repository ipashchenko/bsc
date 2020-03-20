#ifndef BSC__GAIN_H_
#define BSC__GAIN_H_

#include <vector>
#include <set>
#include <Eigen/Core>
#include <valarray>
#include "RNG.h"


using Eigen::MatrixXd;


std::valarray<double> make_normal_random(int number);
std::valarray<double> make_normal_random(int number, double mu, double scale, DNest4::RNG& rng);


class Gain {
    public:
        // Different times for amplitudes and phases of the gains ctor
        Gain(int ant, const std::set<double>& times_amp, const std::set<double>& times_phase);
        // The same times for amplitudes and phases ctor
        Gain(int ant, const std::set<double>& times);

        //Gain(const Gain& other);

        //Gain&operator=(const Gain& other);

        void print(std::ostream &out) const;
        void print_times(std::ostream &out) const;
        void from_prior_amp(DNest4::RNG &rng);
        void from_prior_phase(DNest4::RNG &rng);
        // Generate from prior mean of the GP phase
        void from_prior_phase_mean(DNest4::RNG& rng);
        void from_prior_amp_mean(DNest4::RNG& rng);
        // MH proposals returning logH
        double perturb(DNest4::RNG& rng);
        int size_amp();
        int size_phase();

        const std::valarray<double>& get_amplitudes() const;
        const std::valarray<double>& get_phases() const;
        void sum();

        void print_amplitudes(std::ostream& out) const;
        void print_phases(std::ostream& out) const;
        std::string description() const;


    private:
        // Antenna number
        int ant;
        // Times when gains are measured
        Eigen::VectorXd times_amp;
        Eigen::VectorXd times_phase;
        // Amplitudes and phases of the gains (around means)
        std::valarray<double> amplitudes;
        std::valarray<double> phases;
        // Amplitudes and phases of the gains
        std::valarray<double> amplitudes_full;
        std::valarray<double> phases_full;
        // Mean of the amplitude
        double amp_mean;
        // Mean of the phase
        double phase_mean;

};

#endif //BSC__GAIN_H_
