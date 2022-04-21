//
// Created by admin on 19.04.2022.
//

#include "gtest/gtest.h"
#include "../Kepler.cpp"

TEST(abcd, abcde) {

    const double d =
            static_cast<double>(367 * 2022 - 7 * (2022 + (4 + 9) / 12) / 4 + 275 * 4 / 9 + 12 - 730530) +
            9.0 / 24.0;
    std::cout << d << "\n";

    Position pK = Position(d);
    pK.get_all_Kepler();
    pK.Calculating_Cartesian();
    pK.get_all_Cartesian();

    //Position pD = Position(1.3859e+11, 5.72598e+10, 0, -6.70639e-07, 3.03742e-06, 0); //sun
    //Position pD = Position(1.74319e+10, 4.26557e+10, 1.85987e+09, -5.28306e-06, 2.51763e-06, 5.81689e-07); //mercury
    Position pD = Position(-2.41864e+10, -6.52488e+10, -3.10998e+09, 3.45403e-06, -1.73542e-06, -3.88863e-07);
    pD.get_all_Cartesian();
    pD.Calculating_Kepler();
    pD.get_all_Kepler();
    ASSERT_TRUE(true);
}