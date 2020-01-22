#ifndef BSC__GAIN_H_
#define BSC__GAIN_H_

#include <vector>
#include <set>
#include <Eigen/Core>
#include <valarray>
#include "RNG.h"


using Eigen::MatrixXd;


std::valarray<double> make_normal_random(int number);
std::valarray<double> make_normal_random(int number, DNest4::RNG& rng);


class Gain {
    public:
        // Different times for amplitudes and phases of the gains ctor
        Gain(std::set<double> times_amp, std::set<double> times_phase);
        // The same times for amplitudes and phases ctor
        explicit Gain(std::set<double> times);

        //Gain(const Gain& other);

        //Gain&operator=(const Gain& other);

        void print(std::ostream &out) const;
        void from_prior_global_scale(DNest4::RNG& rng);
        // Generate from prior HP for amplitude and phase
        void from_prior_hp_amp(DNest4::RNG& rng);
        void from_prior_hp_phase(DNest4::RNG& rng);
        // Generate from prior latent variables of amplitude and phase GP
        void from_prior_v_amp(DNest4::RNG &rng);
        void from_prior_v_phase(DNest4::RNG &rng);
        // Generate from prior mean of the GP phase
        void from_prior_phase_mean(DNest4::RNG& rng);
        void from_prior_amp_mean(DNest4::RNG& rng);
        // MH proposals returning logH
        double perturb(DNest4::RNG& rng);
        // Calculate covariance matrixes using current hyperparameters of the gain amplitude and phase
        void calculate_C_amp();
        void calculate_C_phase();
        // Calculate rotation matrixes using current covariance matrixes
        void calculate_L_amp();
        void calculate_L_phase();
        int size_amp();
        int size_phase();
        double get_global_scale() {
            return global_scale;
        }

        // Calculate amplitudes and phases using current values of ``L`` and ``v``.
        void calculate_amplitudes();
        void calculate_phases();
        std::valarray<double>& get_amplitudes();
        std::valarray<double>& get_phases();

        void print_amplitudes(std::ostream& out) const;
        void print_phases(std::ostream& out) const;
        std::string description() const;


    private:
        // Times when gains are measured
        Eigen::VectorXd times_amp;
        Eigen::VectorXd times_phase;
        // Global scale factor
        double global_scale;
        // Hyperparameters of GP for amplitude and phase time dependence
        double logamp_amp;
        double logamp_phase;
        double logscale_amp;
        double logscale_phase;
        // Latent variables ~ N(0, 1)
        std::valarray<double> v_amp;
        std::valarray<double> v_phase;
        // Amplitudes and phases of the gains
        std::valarray<double> amplitudes;
        std::valarray<double> phases;
        // Mean of the amplitude GP
        double amp_mean;
        // Mean of the phase GP
        double phase_mean;
        // Covariance matrix given times and HP
        MatrixXd C_amp;
        MatrixXd C_phase;
        // Choleski-decomposed C that rotates coefficients
        MatrixXd L_amp;
        MatrixXd L_phase;

};

#endif //BSC__GAIN_H_
