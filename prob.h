#include <cassert>
#include <vector>
#include <random>
#include <cmath>

class DirichletMultinomial {
    public:
    DirichletMultinomial(unsigned size, float concentration)
        : K(size), count(size), alpha(concentration) {}

    void Increment(unsigned k) {
        assert(0 <= k && k < K);
        count[k]++;
        N++;
    }

    void Decrement(unsigned k) {
        assert(0 <= k && k < K);
        count[k]--;
        N--;
    }

    float Prob(unsigned k) const {
        assert(0 <= k && k < K);
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
