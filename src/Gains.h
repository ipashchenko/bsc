#ifndef BSC__GAINS_H_
#define BSC__GAINS_H_

#include "Gain.h"
#include "Data.h"
#include "RNG.h"


class Gains {
    public:
        Gains(Data data, int refant);
        //Gains(Gains& other);
        //Gains& operator=(const Gains& other);
        //~Gains();
        //// Setting hyperparameters fo the amplitude and phase GP for all gains
        //void set_hp_phase(std::valarray<double> params);
        //// Setting latent variables for all gains
        //void set_v_amp(std::valarray<double> params);
        //void set_v_phase(std::valarray<double> params);
        // Generate HP for amplitudes and phases for all gains
        void from_prior_hp_amp(DNest4::RNG& rng);
        void from_prior_hp_phase(DNest4::RNG& rng);
        // Generates latent variables from N(0, 1) for all gains
        void from_prior_v_amp();
        void from_prior_v_phase();
        // MH proposal that returns logH
        double perturb(DNest4::RNG& rng);
        int size() const;
        // Getting i-th gain
        Gain& operator[](int i);

        // Calculate covariance matrixes using current hyperparameters for all gains
        void calculate_C_amp();
        void calculate_C_phase();
        // Calculate rotation matrixes using current covariance matrixes for all gains
        void calculate_L_amp();
        void calculate_L_phase();
        // Calculate amplitudes and phases of all gains using current rotation matrixes and latent variables
        void calculate_amplitudes();
        void calculate_phases();
        // Print for all individual antennas gains
        void print_amplitudes(std::ostream& out) const;
        void print_phases(std::ostream& out) const;
        void print_times(std::ostream &out) const;
        void print_hp(std::ostream &out) const;
        void print_v(std::ostream &out) const;
        void print_C(std::ostream &out) const;
        void print_L(std::ostream &out) const;
        std::string description() const;
        void print(std::ostream &out) const;
    private:
        std::vector<Gain> gains;
        int refant{};
        std::vector<int> antennas_changing_gain;
};

#endif //BSC__GAINS_H_
