#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <iostream>
#include <random>
#include <sstream>
#include <iomanip>

namespace serialization {
    
        std::string generate_uuid() {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(0, 15);
            std::uniform_int_distribution<> dis2(8, 11);

            std::stringstream ss;
            ss << std::hex << std::setfill('0');

            // Generate 8-4-4-4-12 UUID format
            for (int i = 0; i < 8; ++i) ss << std::setw(1) << dis(gen);
            ss << "-";
            for (int i = 0; i < 4; ++i) ss << std::setw(1) << dis(gen);
            ss << "-";
            
            // Version 4 UUID marker
            ss << "4"; 
            for (int i = 0; i < 3; ++i) ss << std::setw(1) << dis(gen);
            ss << "-";
            
            // Variant 1 marker (8, 9, a, or b)
            ss << std::setw(1) << dis2(gen); 
            for (int i = 0; i < 3; ++i) ss << std::setw(1) << dis(gen);
            ss << "-";
            
            for (int i = 0; i < 12; ++i) ss << std::setw(1) << dis(gen);

            return ss.str();
    }
};

#endif