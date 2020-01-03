#ifndef BSC__GAINS_H_
#define BSC__GAINS_H_

#include "Gain.h"
#include "Data.h"
#include "RNG.h"


class Gains {
    public:
        //Gains();
        explicit Gains(Data data);
        Gains(Gains& other);
        Gains& operator=(const Gains& other);
        ~Gains();
        // Generate HP for amplitudes and phases for all gains
        void from_prior_hp_amp(DNest4::RNG& rng);
        void from_prior_hp_phase(DNest4::RNG& rng);
        // Generates latent variables from N(0, 1) or use zeros for all gains
        void from_prior_v_amp();
        void from_prior_v_phase();
        // MH proposal that returns logH
        double perturb(DNest4::RNG& rng);
        int size() const;
        // Getting i-th gain
        std::vector<Gain*>& operator[](int i);

        // Calculate covariance matrixes using current hyperparameters for all gains
        void calculate_C_amp();
        void calculate_C_phase();
        // Calculate rotation matrixes using current covariance matrixes for all gains
        void calculate_L_amp();
        void calculate_L_phase();
        // Calculate amplitudes and phases of all gains using current rotation matrixes and latent variables
        void calculate_amplitudes();
        void calculate_phases();
        std::string description() const;
        void print(std::ostream &out) const;
    private:
        // Vector length ``n_ant`` of vectors length ``n_IF`` each
        std::vector<std::vector<Gain*>> gains;
};

#endif //BSC__GAINS_H_
