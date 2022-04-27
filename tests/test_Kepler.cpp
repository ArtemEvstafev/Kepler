//
// Created by admin on 19.04.2022.
//

#include "gtest/gtest.h"
#include "../Kepler.cpp"

TEST(KeplerTests, Mercury) {

    std::cout<<"FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position MercuryK = Position();
    MercuryK.set_kepler(48.3313, 7.0047, 29.1241, 0.387098, 0.205635, 168.6562);
    MercuryK.get_all_kepler();
    MercuryK.get_all_cartesian();

    std::cout<<"FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position MercuryC = Position();
    MercuryC.set_cartesian( -2.41873e+10, -6.52484e+10, -3.10986e+09, 35886.8, -14519, -4479.77);
    MercuryC.get_all_cartesian();
    MercuryC.get_all_kepler();

    ASSERT_TRUE(true);
}
TEST(KeplerTests, Sun){

    std::cout<<"FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position SunK = Position();
    SunK.set_kepler(0.0, 0.0, 282.9404, 1.000000, 0.016709, 356.0470);
    SunK.get_all_kepler();
    SunK.get_all_cartesian();

    std::cout<<"FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position SunC = Position();
    SunC.set_cartesian(2.26383e+10, -1.45352e+11, 0.0, 29919.1, 4695.76, 0.0);
    SunC.get_all_cartesian();
    SunC.get_all_kepler();

    ASSERT_TRUE(true);
}

TEST(KeplerTests, Mars){

    std::cout<<"FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position MarsK = Position();
    MarsK.set_kepler(49.5574, 1.8497, 286.5016, 1.523688, 0.093405, 18.6021);
    MarsK.get_all_kepler();
    MarsK.get_all_cartesian();

    std::cout<<"FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position MarsC = Position();
    MarsC.set_cartesian(2.07859e+11, -5.40842e+09, -5.22205e+09, 1561.45, 26290.5, 512.383);
    MarsC.get_all_cartesian();
    MarsC.get_all_kepler();

    ASSERT_TRUE(true);
}