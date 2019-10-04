#ifndef BSC__GAIN_H_
#define BSC__GAIN_H_

#include <vector>
#include <set>
#include <Eigen/Core>
#include <valarray>
#include "RNG.h"


using Eigen::MatrixXd;


std::valarray<double> make_normal_random(int number);


class BasicGain {
    public:
        virtual BasicGain* clone() const = 0;
        virtual void from_prior(DNest4::RNG& rng) = 0;
        virtual double perturb(DNest4::RNG& rng) = 0;
        std::valarray<double>& get_amplitudes() {
            return amplitudes;
        }
        std::valarray<double>& get_phases() {
            return phases;
        }
        //void print_amplitudes(std::ostream& out) const;
        //void print_phases(std::ostream& out) const;
        virtual std::string description() const = 0;
        virtual void print(std::ostream &out) const = 0;
    protected:
        // Amplitudes and phases of the gains
        std::valarray<double> amplitudes;
        std::valarray<double> phases;
        // Times when gains are measured
        Eigen::VectorXd times_amp;
        Eigen::VectorXd times_phase;
};


// https://stackoverflow.com/a/5027623
template<class Derived>
class Cloneable : public BasicGain {
    virtual BasicGain* clone() const {
        return new Derived(static_cast<const Derived&>(*this));
    }
};


class Gain : public Cloneable<Gain> {
    public:
        // Different times for amplitudes and phases of the gains ctor
        Gain(std::set<double> times_amp, std::set<double> times_phase);
        // The same times for amplitudes and phases ctor
        explicit Gain(std::set<double> times);

        Gain(const Gain& other);

        //Gain&operator=(const Gain& other);

        std::string description() const override ;
        void print(std::ostream &out) const override ;
        //void print_times(std::ostream &out) const;
        //void print_hp(std::ostream &out) const;
        //void print_v(std::ostream &out) const;
        //void print_C(std::ostream &out) const;
        //void print_L(std::ostream &out) const;
        // Set hyperparameters of amplitude and phase GP.
        //void set_hp_amp(std::valarray<double> params);
        //void set_hp_phase(std::valarray<double> params);
        // Set latent variables of amplitude and phase GP
        //void set_v_amp(std::valarray<double> params);
        //void set_v_phase(std::valarray<double> params);
        void from_prior(DNest4::RNG& rng) override ;
        // Generate from prior HP for amplitude and phase
        void from_prior_hp_amp(DNest4::RNG& rng);
        void from_prior_hp_phase(DNest4::RNG& rng);
        // Generate from prior latent variables of amplitude and phase GP
        void from_prior_v_amp();
        void from_prior_v_phase();
        // Generate from prior mean of the GP phase
        void from_prior_phase_mean(DNest4::RNG& rng);
        // MH proposals returning logH
        double perturb(DNest4::RNG& rng) override ;
        // Calculate covariance matrixes using current hyperparameters of the gain amplitude and phase
        void calculate_C_amp();
        void calculate_C_phase();
        // Calculate rotation matrixes using current covariance matrixes
        void calculate_L_amp();
        void calculate_L_phase();
        int size_amp();
        int size_phase();

        // Calculate amplitudes and phases using current values of ``L`` and ``v``.
        void calculate_amplitudes();
        void calculate_phases();


    private:
        // Hyperparameters of GP for amplitude and phase time dependence
        double logamp_amp;
        double logamp_phase;
        double logscale_amp;
        double logscale_phase;
        // Latent variables ~ N(0, 1)
        std::valarray<double> v_amp;
        std::valarray<double> v_phase;
        // Mean of the phase GP
        double phase_mean;
        // Covariance matrix given times and HP
        MatrixXd C_amp;
        MatrixXd C_phase;
        // Choleski-decomposed C that rotates coefficients
        MatrixXd L_amp;
        MatrixXd L_phase;
};


class ConstantGain : public Cloneable<ConstantGain> {
    public:

        // Different times for amplitudes and phases of the gains ctor
        ConstantGain(std::set<double> times_amp, std::set<double> times_phase);
        // The same times for amplitudes and phases ctor
        explicit ConstantGain(std::set<double> times);
        ConstantGain(const ConstantGain& other);

        void from_prior(DNest4::RNG& rng) override ;
        double perturb(DNest4::RNG& rng) override ;
        std::string description() const override ;
        void print(std::ostream &out) const override ;
    private:
        // This two always return arrays of constant values
        // Constant over all observing time amplitude and phase of the gain
        double amplitude;
        double phase;
        // Just templates to be modified by one constant value
        std::valarray<double> amplitudes_template;
        std::valarray<double> phases_template;
        // Arrays of repeated values for consistency
        std::valarray<double> amplitudes;
        std::valarray<double> phases;
};


#endif //BSC__GAIN_H_
