//
// Created by artem on 18.04.22.
//
#include "Kepler.hpp"
#include "overloads.hpp"

//template<typename T>
class Position {
private:
    double N; //Longitude of the ascending node
    double i; //Inclination
    double w; //Argument of periapsis
    double a; //Semimajor axis
    double e; //Eccentricity
    double M; //Mean anomaly

    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    //constatnts
    double m = 1.32712440019E+20;
    double AU = 149597870700.0;
public:

    Position() {};

    ~Position() {};

    void set_cartesian(double x, double y, double z, double vx, double vy, double vz) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->vx = vx;
        this->vy = vy;
        this->vz = vz;
        calculating_kepler();
    }

    void set_kepler(double n, double i, double w, double a, double e, double m) {
        this->N = rad(n);
        this->i = rad(i);
        this->w = rad(w);
        this->a = a * AU;
        this->e = e;
        this->M = rad(m);
        calculating_cartesian();
    }


    double rad(double deg) {
        return deg / 180.0 * M_PI;
    }

    double calcE(const double M, double e) {
        double E = M;
        while (std::abs(E - M - e * sin(E)) > 0.00001)
            E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
        return E;
    }

    void calculating_cartesian() {

        double E = calcE(M, e);
        double v = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
        double r = a * (1 - e * cos(E));

        double ox = r * cos(v);
        double oy = r * sin(v);
        double oz = 0;
        double ovx = -sqrt(m * a) / r * sin(E);
        double ovy = sqrt(m * a) / r * sqrt(1 - e * e) * cos(E);
        double ovz = 0;

        x = ox * (cos(w) * cos(N) - sin(w) * cos(i) * sin(N)) - oy * (sin(w) * cos(N) + cos(w) * cos(i) * sin(N));
        y = ox * (cos(w) * sin(N) + sin(w) * cos(i) * cos(N)) + oy * (cos(w) * cos(i) * cos(N) - sin(w) * sin(N));
        z = ox * (sin(w) * sin(i)) + oy * (cos(w) * sin(i));

        vx = ovx * (cos(w) * cos(N) - sin(w) * cos(i) * sin(N)) - ovy * (sin(w) * cos(N) + cos(w) * cos(i) * sin(N));
        vy = ovx * (cos(w) * sin(N) + sin(w) * cos(i) * cos(N)) + ovy * (cos(w) * cos(i) * cos(N) - sin(w) * sin(N));
        vz = ovx * (sin(w) * sin(i)) + ovy * (cos(w) * sin(i));

    };

    template<typename T>
    std::vector<T> vectmul(const std::vector<T> &c, const std::vector<T> &b) {
        std::vector<T> result(c.size());
        result[0] = c[1] * b[2] - b[1] * c[2];
        result[1] = -c[0] * b[2] + b[0] * c[2];
        result[2] = c[0] * b[1] - b[0] * c[1];
        return result;
    }


    void calculating_kepler() {

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

        if (e != 1)
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
            M += rad(360);

        a = 1 / (2 / norm(r) - norm(dr) * norm(dr) / m);

    };

    void get_all_kepler() {
        std::cout << "KEPLER ELEMENTS\n" << "N " << N << "\n" << "i " << i << "\n" << "w " << w << "\n" << "a " << a
                  << "\n" << "e " << e << "\n" << "M " << M << "\n" << "(" << N << ", " << i << ", " << w
                  << ", " << a << ", " << e << ", " << M << ")\n";;
    };

    void get_all_cartesian() {
        std::cout << "CARTESIAN COORDINATES\n" << "x " << x << "\n" << "y " << y << "\n" << "z " << z << "\n" << "vx "
                  << vx << "\n" << "vy " << vy << "\n" "vz " << vz << "\n" << "(" << x << ", " << y << ", " << z
                  << ", " << vx << ", " << vy << ", " << vz << ")\n";
    };

};


