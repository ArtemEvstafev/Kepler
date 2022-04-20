//
// Created by admin on 19.04.2022.
//

#include "gtest/gtest.h"
#include "../Kepler.cpp"

TEST(abcd, abcde) {

    const double d =
            static_cast<double>(367 * 2022 - 7 * (2022 + (4 + 9) / 12) / 4 + 275 * 4 / 9 + 12 - 730530) +
            9.0 / 24.0;
    const double ecl = 23.4393 - 3.563E-7 * d;
    std::cout << d << "\t" << ecl << "\n";
    Position p = Position(d, ecl);
    p.get_all_Kepler();
    p.Calculating_Cartesian();
    p.get_all_Cartesian();
    ASSERT_TRUE(true);
}