#include <cassert>
#include <vector>
#include <random>
#include <cmath>
#include <mutex>
#include <iostream>

/* A thread-safe Multinomial distribution with a Dirichlet prior */

struct DirichletMultinomial {
    DirichletMultinomial(unsigned size, float concentration)
        : K(size), count(size), N(0), alpha(concentration) {}

    void Increment(unsigned k) {
        assert(k < K);
        std::lock_guard<std::mutex> guard(count_lock);
        count[k]++;
        N++;
    }

    void Decrement(unsigned k) {
        assert(k < K);
        std::lock_guard<std::mutex> guard(count_lock);
        count[k]--;
        N--;
    }

    float Prob(unsigned k) const { // Posterior predictive: p(x_n=k | x^-n)
        assert(k < K);
        return (alpha + count[k])/(K * alpha + N);
    }

    double LogLikelihood() const { // p(x|alpha) = \int_theta p(x|theta) p(theta|alpha)
        double ll = lgamma(K * alpha) - K * lgamma(alpha) - lgamma(K * alpha + N);
        for(unsigned k = 0; k < K; k++)
            ll += lgamma(alpha + count[k]);
        return ll;
    }

    unsigned K;
    float alpha;
    unsigned N;
    std::vector<unsigned> count;
    std::mutex count_lock;

    friend std::ostream& operator<<(std::ostream&, const DirichletMultinomial&);
};


std::ostream& operator<<(std::ostream& os, const DirichletMultinomial& m) {
    unsigned support = 0;
    for(unsigned k = 0; k < m.K; k++)
        support += (m.count[k] > 0);
    return os << "Multinomial(N=" << m.N << " |support|=" << support << ")"
        << " ~ Dir(K=" << m.K << ", alpha=" << m.alpha << ")";
}

/* A Geometric distribution with a Beta prior */

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

    float Stop() const { // E[p|alpha] - used to approximate posterior predictive
        return (alpha + N)/(alpha + N + beta + L); // mean = 1/p - 1 = (beta+L)/(alpha+N)
    }

    float Prob(unsigned l) const {
        float p = (alpha + N) / (alpha + N + beta + L);
        for(int k = 0; k <= l-1; k++)
            p *= (beta + L + l)/(alpha + N + 1 + beta + L + l);
        return p;
    }

    double LogLikelihood() const {
        return lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) 
            + lgamma(alpha + N) + lgamma(beta + L) - lgamma(alpha + N + beta + L);
    }

    friend std::ostream& operator<<(std::ostream&, const BetaGeometric&);
};


std::ostream& operator<<(std::ostream& os, const BetaGeometric& m) {
    return os << "Geometric(N=" << m.N << ", L=" << m.L << ")"
        <<" ~ Beta(" << m.alpha << ", " << m.beta << ")";
}


/* Utility functions for generating random numbers */

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
