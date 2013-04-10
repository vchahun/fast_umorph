#include <fst/fstlib.h>
#include "vocabulary.h"
#include "corpus.h"
#include "prob.h"
#include "trie.h"
#include "banana.h"
#include "pss_model.h"

const std::string FormatSegmentation(const Segmentation& seg,
        const Vocabulary& substring_vocabulary) {
    std::string res;
    for(unsigned p: seg.prefixes)
        res += substring_vocabulary.Convert(p) + "^";
    res += "<" + substring_vocabulary.Convert(seg.stem) + ">";
    for(unsigned s: seg.suffixes)
        res += "^" + substring_vocabulary.Convert(s);
    return res;
}

int main(int argc, char** argv) {
    assert(argc == 2);

    const float alpha_prefix = 0.001;
    const float alpha_stem = 0.001;
    const float alpha_suffix = 0.001;

    Vocabulary word_vocabulary;
    Vocabulary substring_vocabulary;

    Corpus corpus(std::cin, word_vocabulary);
    std::cerr << "Read " << corpus.Size() << " tokens, "
        << word_vocabulary.Size() << " types\n";


    std::vector<Trie> tries;
    for(const std::string& word: word_vocabulary) {
        Trie trie;
        for(unsigned i = 0; i <= word.size(); i++) {
            for(unsigned j = 1; i + j <= word.size(); j++) {
                const std::string substring = word.substr(i, j);
                trie.Insert(substring, substring_vocabulary.Encode(substring));
            }
        }
        tries.push_back(trie);
    }

    std::cerr << "Found " << substring_vocabulary.Size() << " substrings\n";

    SegmentationModel model(alpha_prefix, alpha_stem, alpha_suffix,
           word_vocabulary, substring_vocabulary.Size(), tries);

    std::random_device rd;
    std::mt19937 engine(rd());

    const unsigned n_iterations = atoi(argv[1]);

    std::vector<Segmentation> segs;
    for(unsigned it = 0; it < n_iterations; it++) {
        unsigned wid = 0;
        for(auto& sentence: corpus) {
            for(auto& word: sentence) {
                if(it > 0) model.Decrement(word, segs[wid]);
                const Segmentation seg = model.Increment(word, engine, (it==0));
                if(it > 0) segs[wid] = seg;
                else segs.push_back(seg);
                wid++;
                if(wid % 100 == 0) std::cerr << ".";
            }
        }
        std::cerr << "\n";
        std::cerr << model << "\n";
        if(it % 10 == 9) {
            std::cerr << "Iteration " << (it+1) << "/" << n_iterations << "\n";
            double ll = model.LogLikelihood();
            double ppl = exp(-ll/corpus.Tokens());
            std::cerr << "LL=" << ll << " ppl=" << ppl << "\n";
        }
    }

    for(unsigned w = 0; w < word_vocabulary.Size(); w++) {
        const std::string& word = word_vocabulary.Convert(w);
        const Segmentation seg = model.Decode(w);
        std::cout << word << "\t" << FormatSegmentation(seg, substring_vocabulary) << "\n";
    }
}
