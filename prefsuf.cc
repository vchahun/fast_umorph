#include "prob.h"
#include "vocabulary.h"
#include "corpus.h"

class LexiconModel {
    public:
    LexiconModel(float alpha_t, float alpha_f,
            const Vocabulary& word_vocabulary,
            const Vocabulary& prefix_vocabulary,
            const Vocabulary& suffix_vocabulary) :
        word_vocabulary(word_vocabulary),
        prefix_vocabulary(prefix_vocabulary),
        suffix_vocabulary(suffix_vocabulary),
        prefix_model(prefix_vocabulary.Size(), alpha_t),
        suffix_model(suffix_vocabulary.Size(), alpha_f) {}

    unsigned Increment(unsigned w, std::mt19937& engine, bool initialize=false) {
        const std::string& word = word_vocabulary.Convert(w);
        unsigned t, f, split;
        if(initialize) {
            split = prob::randint(engine, 0, word.size());
            t = prefix_vocabulary.Convert(word.substr(0, split));
            f = suffix_vocabulary.Convert(word.substr(split));
        }
        else {
            float x = prob::random(engine) * Prob(w);
            float analysis_prob;
            for(split = 0 ; split <= word.size(); split++) {
                t = prefix_vocabulary.Convert(word.substr(0, split));
                f = suffix_vocabulary.Convert(word.substr(split));
                analysis_prob = prefix_model.Prob(t) * suffix_model.Prob(f);
                if(x < analysis_prob || split == word.size()) break;
                x -= analysis_prob;
            }
        }
        prefix_model.Increment(t);
        suffix_model.Increment(f);
        return split;
    }

    void Decrement(unsigned w, unsigned split) {
        const std::string& word = word_vocabulary.Convert(w);
        unsigned t = prefix_vocabulary.Convert(word.substr(0, split));
        unsigned f = suffix_vocabulary.Convert(word.substr(split));
        prefix_model.Decrement(t);
        suffix_model.Decrement(f);
    }
    
    float Prob(unsigned w) const {
        const std::string& word = word_vocabulary.Convert(w);
        float prob = 0;
        for(unsigned split = 0; split <= word.size(); split++) {
            unsigned t = prefix_vocabulary.Convert(word.substr(0, split));
            unsigned f = suffix_vocabulary.Convert(word.substr(split));
            prob += prefix_model.Prob(t) * suffix_model.Prob(f);
        }
        return prob;
    }

    unsigned Decode(unsigned w) const {
        const std::string& word = word_vocabulary.Convert(w);
        float max_prob = -1;
        unsigned best_split;
        for(unsigned split = 0; split <= word.size(); split++) {
            unsigned t = prefix_vocabulary.Convert(word.substr(0, split));
            unsigned f = suffix_vocabulary.Convert(word.substr(split));
            float prob = prefix_model.Prob(t) * suffix_model.Prob(f);
            if(prob >= max_prob) {
                max_prob = prob;
                best_split = split;
            }
        }
        return best_split;
    
    }

    double LogLikelihood() const {
        return prefix_model.LogLikelihood() + suffix_model.LogLikelihood();
    }

    private:
    const Vocabulary& word_vocabulary, prefix_vocabulary, suffix_vocabulary;
    DirichletMultinomial prefix_model, suffix_model;
};

int main(int argc, char** argv) {
    assert(argc == 2);

    Vocabulary word_vocabulary, prefix_vocabulary, suffix_vocabulary;
    Corpus corpus(std::cin, word_vocabulary);
    std::cerr << "Read " << corpus.Size() << " tokens, "
        << word_vocabulary.Size() << " types\n";

    for(auto& seg: corpus) {
        for(auto& w: seg) {
            const std::string& word = word_vocabulary.Convert(w);
            unsigned t, f;
            for(unsigned split = 0; split <= word.size(); split++) {
                prefix_vocabulary.Encode(word.substr(0, split));
                suffix_vocabulary.Encode(word.substr(split));
            }
        }
    }

    std::cerr << "Found " << prefix_vocabulary.Size() << " prefixes, "
        << suffix_vocabulary.Size() << " suffixes\n";

    LexiconModel model(0.001, 0.001, word_vocabulary, prefix_vocabulary, suffix_vocabulary);
    std::random_device rd;
    std::mt19937 engine(rd());

    const unsigned n_iterations = atoi(argv[1]);

    std::vector<unsigned> splits;
    for(unsigned it = 0; it < n_iterations; it++) {
        unsigned wid = 0;
        for(auto& segment: corpus) {
            for(auto& word: segment) {
                if(it > 0) model.Decrement(word, splits[wid]);
                unsigned split = model.Increment(word, engine, (it==0));
                if(it > 0) splits[wid] = split;
                else splits.push_back(split);
                wid++;
            }
        }
        if(it % 10 == 9) {
            std::cerr << "Iteration " << (it+1) << "/" << n_iterations << "\n";
            double ll = model.LogLikelihood();
            double ppl = exp(-ll/corpus.Tokens());
            std::cerr << "LL=" << ll << " ppl=" << ppl << "\n";
        }
    }

    for(unsigned w = 0; w < word_vocabulary.Size(); w++) {
        const std::string& word = word_vocabulary.Convert(w);
        unsigned split = model.Decode(w);
        std::cout << word << "\t_\t"
            << word.substr(0, split) << "\t"
            << word.substr(split) << "\n";
    }
}
