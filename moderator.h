#pragma once
#ifndef CHERNOBYL_CORE_MODERATOR_H
#define CHERNOBYL_CORE_MODERATOR_H

#include "particle.h"

class Moderator : public Particle {
public:
    std::vector<FissionSplit> bohr_wheeler(std::optional<float> energy_input = std::nullopt) const {
        return {
            FissionSplit{ mass, 0, 0, 1.0f }
        };
    }
};
#endif //CHERNOBYL_CORE_MODERATOR_H