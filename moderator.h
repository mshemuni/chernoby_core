//
// Created by niaei on 11.01.2026.
//

#pragma once
#ifndef CHERNOBY_CORE_MODERATOR_H
#define CHERNOBY_CORE_MODERATOR_H

#include "particle.h"

class Moderator : public Particle {
public:
    std::vector<FissionSplit> bohr_wheeler(std::optional<float> energy_input = std::nullopt) const {
        return {
            FissionSplit{ mass, 0, 0, 1.0f }
        };
    }
};
#endif //CHERNOBY_CORE_MODERATOR_H