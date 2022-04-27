//
// Created by admin on 20.04.2022.
//
template<typename T>
T operator*(const std::vector<T> &a, const std::vector<T> &b) {
    T result = 0;
    for (int i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<T> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }

    return result;
}

template<typename T>
std::vector<T> operator*(const T &alpha, const std::vector<T> &b) {
    std::vector<T> result(b.size());
    for (int i = 0; i < result.size(); i++) {
        result[i] = alpha * b[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator*(const std::vector<T> &b, const T &alpha) {
    std::vector<T> result(b.size());
    for (int i = 0; i < result.size(); i++) {
        result[i] = alpha * b[i];
    }
    return result;
}

template<typename T>
T norm(const std::vector<T> &vector) {
    T norm = static_cast<T>(0);
    for (const auto &elm: vector) {
        norm += elm * elm;
    }
    norm = std::sqrt(norm);
    return norm;
}
template<typename T>
std::ostream& operator<<(std::ostream &os, std::vector<T>& vec) {
    std::cout<<"(";
    for(auto& e: vec)
        std::cout<<e<<" ";
    std::cout<<")\n";
    return os;
}

