#include <cassert>
#include <vector>
#include <random>
#include <cmath>

struct DirichletMultinomial {
    DirichletMultinomial(unsigned size, float concentration)
        : K(size), count(size), alpha(concentration) {}

    void Increment(unsigned k) {
        assert(k < K);
        count[k]++;
        N++;
    }

    void Decrement(unsigned k) {
        assert(k < K);
        count[k]--;
        N--;
    }

    float Prob(unsigned k) const {
        assert(k < K);
        return (alpha + count[k])/(K * alpha + N);
    }

    double LogLikelihood() const {
        double ll = lgamma(K * alpha) - lgamma(K * alpha + N) - K * lgamma(alpha);
        for(unsigned k = 0; k < K; k++)
            ll += lgamma(alpha + count[k]);
        return ll;
    }

    unsigned K;
    float alpha;
    unsigned N;
    std::vector<unsigned> count;

    friend std::ostream& operator<<(std::ostream&, const DirichletMultinomial&);
};


std::ostream& operator<<(std::ostream& os, const DirichletMultinomial& m) {
    unsigned support = 0;
    for(unsigned k = 0; k < m.K; k++)
        support += (m.count[k] > 0);
    return os << "Multinomial(N=" << m.N << " |support|=" << support << ")"
        << " ~ Dir(K=" << m.K << ", alpha=" << m.alpha << ")";
}

struct BetaGeometric {
    unsigned L, N;
    float alpha, beta;

    BetaGeometric(float alpha, float beta) : L(0), N(0), alpha(alpha), beta(beta) {}

    void Increment(unsigned l) {
        L += l;
        N++;
    }

    void Decrement(unsigned l){
        L -= l;
        N--;
    }

    float Stop() const {
        return (alpha + N)/(alpha + beta + N + L);
    }

    float Prob(unsigned l) const {
        float p = Stop();
        return p * pow(1 - p, l); // mean = 1/p - 1
    }

    double LogLikelihood() const {
        return -1;
    }

    friend std::ostream& operator<<(std::ostream&, const BetaGeometric&);
};


std::ostream& operator<<(std::ostream& os, const BetaGeometric& m) {
    return os << "Geometric(N=" << m.N << ", L=" << m.L << ")"
        <<" ~ Beta(" << m.alpha << ", " << m.beta << ")";
}


namespace prob {

template<typename Engine> 
double random(Engine& engine) {
    return std::uniform_real_distribution<>(0, 1)(engine);
}

template<typename Engine> 
double randint(Engine& engine, int a, int b) {
    return std::uniform_int_distribution<>(a, b)(engine);
}

}
