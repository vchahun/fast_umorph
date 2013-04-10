#include <cassert>
#include <vector>
#include <random>
#include <cmath>

class DirichletMultinomial {
    public:
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

    private:
    unsigned K;
    float alpha;
    unsigned N;
    std::vector<unsigned> count;
};

struct BetaGeometric {
    unsigned L, N;
    float alpha, beta;

    BetaGeometric(float alpha, float beta) : L(0), N(0), alpha(alpha), beta(beta) {}

    void Increment(unsigned l) {
        L += l;
        N += 1;
    }

    void Decrement(unsigned l ){
        L -= l;
        N += 1;
    }

    float Stop() const {
        return (alpha + N)/(alpha + beta + N + L);
    }

    float Prob(unsigned l) const {
        float p = Stop();
        return p * pow(1 - p, l);
    }

    double LogLikelihood() const {
        return -1;
    }
};


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
