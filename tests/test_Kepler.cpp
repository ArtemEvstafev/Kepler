//
// Created by admin on 19.04.2022.
//

#include "gtest/gtest.h"
#include "../Kepler.cpp"

const double PERCENT = 0.01;    //точность в процентах

TEST(KeplerTests, Mercury) {

    std::vector<double> res;
    std::vector<double> real;
    std::cout << "FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position MercuryK = Position();
    MercuryK.set_kepler(0.84354, 0.122255, 0.508311, 0.387098, 0.205635, 2.94361);

    real = {-2.41873e+10, -6.52484e+10, -3.10986e+09, 35886.8, -14519, -4479.77};
    res = MercuryK.get_all_cartesian();

    for (int i = 0; i < real.size(); i++) {
        ASSERT_TRUE(std::abs(res[i] - real[i]) <= std::abs(real[i] * PERCENT));
    }

    std::cout << "FROM CARTESIAN STATE VECTORS TO KEPLER ELEMENTS\n";
    Position MercuryC = Position();
    MercuryC.set_cartesian(-2.41873e+10, -6.52484e+10, -3.10986e+09, 35886.8, -14519, -4479.77);

    real = {0.84354, 0.122255, 0.508311, 0.387098, 0.205635, 2.94361};
    res = MercuryC.get_all_kepler();

    for (int i = 0; i < real.size(); i++) {
        ASSERT_TRUE(std::abs(res[i] - real[i]) <= std::abs(real[i] * PERCENT));
    }
}

TEST(KeplerTests, Sun) {

    /*
     * Тест с Солнцем не проходит, так как углы вырождаются и нельзя однозначно определить параметры
     */
    std::vector<double> res;
    std::vector<double> real;
    std::cout << "FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position SunK = Position();
    SunK.set_kepler(0.0, 0.0, 4.93824, 1.000000, 0.016709, 6.21419);

    real = {2.26383e+10, -1.45352e+11, 0.0, 29919.1, 4695.76, 0.0};
    res = SunK.get_all_cartesian();
    for (int i = 0; i < real.size(); i++) {
        ASSERT_TRUE(std::abs(res[i] - real[i]) <= std::abs(real[i] * PERCENT));
    }

    std::cout << "FROM CARTESIAN STATE VECTORS TO KEPLER ELEMENTS\n";
    Position SunC = Position();
    SunC.set_cartesian(2.26383e+10, -1.45352e+11, 0.0, 29919.1, 4695.76, 0.0);

    real = {0.0, 0.0, 4.93824, 1.000000, 0.016709, 6.21419};
    res = SunC.get_all_kepler();
    std::cout << res;
    for (int i = 0; i < real.size(); i++) {
        std::cout << std::abs(res[i] - real[i]) << " <= " << std::abs(real[i] * PERCENT) << '\n';
        ASSERT_TRUE(std::abs(res[i] - real[i]) <= std::abs(real[i] * PERCENT));
    }
}

TEST(KeplerTests, Mars) {

    std::vector<double> res;
    std::vector<double> real;
    std::cout << "FROM KEPLER ELEMENTS TO CARTESIAN STATE VECTORS\n";
    Position MarsK = Position();
    MarsK.set_kepler(0.86494, 0.0322834, 5.0004, 1.523688, 0.093405, 0.324668);

    real = {2.07859e+11, -5.40842e+09, -5.22205e+09, 1561.45, 26290.5, 512.383};
    res = MarsK.get_all_cartesian();

    for (int i = 0; i < real.size(); i++) {
        ASSERT_TRUE(std::abs(res[i] - real[i]) <= std::abs(real[i] * PERCENT));
    }

    std::cout << "FROM CARTESIAN STATE VECTORS TO KEPLER ELEMENTS\n";
    Position MarsC = Position();
    MarsC.set_cartesian(2.07859e+11, -5.40842e+09, -5.22205e+09, 1561.45, 26290.5, 512.383);

    real = {0.86494, 0.0322834, 5.0004, 1.523688, 0.093405, 0.324668};
    res = MarsC.get_all_kepler();

    for (int i = 0; i < real.size(); i++) {
        ASSERT_TRUE(std::abs(res[i] - real[i]) <= std::abs(real[i] * PERCENT));
    }
}