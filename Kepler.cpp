//
// Created by artem on 18.04.22.
//
#include "Kepler.hpp"
#include "overloads.hpp"


Position::Position() {};

Position::~Position() {};

void Position::set_cartesian(double x, double y, double z, double vx, double vy, double vz) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->vx = vx;
    this->vy = vy;
    this->vz = vz;
    calculating_kepler();
}

void Position::set_kepler(double n, double i, double w, double a, double e, double m) {
    this->N = n;
    this->i = i;
    this->w = w;
    this->a = a * AU;
    this->e = e;
    this->M = m;
    calculating_cartesian();
}

double Position::calcE(const double M, const double e) {
    double E = M;
    while (std::abs(E - M - e * sin(E)) > tolerance)
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
    return E;
}

void Position::calculating_cartesian() {

    double E = calcE(M, e);
    double v = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
    double r = a * (1 - e * cos(E));

    double ox = r * cos(v);
    double oy = r * sin(v);
    double oz = 0;
    double ovx = -sqrt(m * a) / r * sin(E);
    double ovy = sqrt(m * a) / r * sqrt(1 - e * e) * cos(E);
    double ovz = 0;
    double cosw = cos(w);
    double cosN = cos(N);
    double cosi = cos(i);
    double sinw = sin(w);
    double sinN = sin(N);
    double sini = sin(i);
    x = ox * (cosw * cosN - sinw * cosi * sinN) - oy * (sinw * cosN + cosw * cosi * sinN);
    y = ox * (cosw * sinN + sinw * cosi * cosN) + oy * (cosw * cosi * cosN - sinw * sinN);
    z = ox * (sinw * sini) + oy * (cosw * sini);

    vx = ovx * (cosw * cosN - sinw * cosi * sinN) - ovy * (sinw * cosN + cosw * cosi * sinN);
    vy = ovx * (cosw * sinN + sinw * cosi * cosN) + ovy * (cosw * cosi * cosN - sinw * sinN);
    vz = ovx * (sinw * sini) + ovy * (cosw * sini);

};


void Position::calculating_kepler() {

    std::vector<double> r = {x, y, z};
    std::vector<double> dr = {vx, vy, vz};
    std::vector<double> h = vectmul(r, dr);
    std::vector<double> ve = 1 / m * vectmul(dr, h) - 1 / norm(r) * r;
    std::vector<double> n = {-h[1], h[0], 0.0};
    double v;

    if (r * dr >= 0) {
        v = acos(ve * r / norm(ve) / norm(r));
    } else
        v = 2 * M_PI - acos(ve * r / norm(ve) / norm(r));

    if (norm(h) != 0) {
        i = acos(h[2] / norm(h));
    } else
        i = 0.0;
    e = norm(ve);
    double E;

    if (std::abs(e - 1.0) > tolerance)
        //if (e != 1)
        E = 2 * atan(tan(v / 2) / sqrt((1 + e) / (1 - e)));
    else
        E = 0.0;

    if (n[1] >= 0) {
        if (norm(n) > 0)
            N = acos(n[0] / norm(n));
        else
            N = 0.0;
    } else
        N = 2 * M_PI - acos(n[0] / norm(n));

    if (ve[2] >= 0) {
        if (norm(n) > 0 && norm(ve) > 0)
            w = acos(ve * n / norm(ve) / norm(n));
        else
            w = 0.0;
    } else
        w = 2 * M_PI - acos(ve * n / norm(ve) / norm(n));

    M = E - e * sin(E);
    if (M < 0)
        M += 2 * M_PI;

    a = 1 / (2 / norm(r) - norm(dr) * norm(dr) / m);

};

std::vector<double> Position::get_all_kepler() {
    /*std::cout << "KEPLER ELEMENTS\n" << "N " << N << "\n" << "i " << i << "\n" << "w " << w << "\n" << "a " << a
              << "\n" << "e " << e << "\n" << "M " << M << "\n" << "(" << N << ", " << i << ", " << w
              << ", " << a << ", " << e << ", " << M << ")\n";*/
    std::vector<double> res = {N, i, w, a / AU, e, M};
    return res;
};

std::vector<double> Position::get_all_cartesian() {
    /*std::cout << "CARTESIAN COORDINATES\n" << "x " << x << "\n" << "y " << y << "\n" << "z " << z << "\n" << "vx "
              << vx << "\n" << "vy " << vy << "\n" "vz " << vz << "\n" << "(" << x << ", " << y << ", " << z
              << ", " << vx << ", " << vy << ", " << vz << ")\n";*/
    std::vector<double> res = {x, y, z, vx, vy, vz};
    return res;
};



