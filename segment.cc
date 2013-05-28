#include <fst/fstlib.h>
#include <thread>
#include "thread_pool.h"
#include "vocabulary.h"
#include "corpus.h"
#include "prob.h"
#include "trie.h"
#include "banana.h"
#include "pss_model.h"

const unsigned NTHREADS = 8;

const std::string FormatSegmentation(const Segmentation& seg,
        const Vocabulary& substring_vocabulary,
        const string morpheme_separator = "^",
        const string prefix_separator = "<",
        const string suffix_separator = ">") {
    std::string res;
    for(unsigned p: seg.prefixes)
        res += substring_vocabulary.Convert(p) + morpheme_separator;
    res += prefix_separator + substring_vocabulary.Convert(seg.stem) + suffix_separator;
    for(unsigned s: seg.suffixes)
        res += morpheme_separator + substring_vocabulary.Convert(s);
    return res;
}

int main(int argc, char** argv) {
    if(argc != 5) {
        std::cerr << "Usage: "
            << argv[0] << " n_iter alpha_prefix alpha_stem alpha_suffix\n";
        exit(1);
    }

    const unsigned n_iterations = atoi(argv[1]);
    const float alpha_prefix = atof(argv[2]);
    const float alpha_stem = atof(argv[3]);
    const float alpha_suffix = atof(argv[4]);

    Vocabulary word_vocabulary;
    Vocabulary substring_vocabulary;

    /* Read vocabulary from standard input */
    Corpus corpus(std::cin, word_vocabulary);
    std::cerr << "Read " << corpus.Size() << " sentences, "
        << corpus.Tokens() << " tokens, "
        << word_vocabulary.Size() << " types\n";

    /* Create substring tries which are used as a based to build
     * segmentation lattices */
    std::vector<Trie> tries;
    for(const std::string& word: word_vocabulary) {
        CheckChars(word);
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

    /* Initialize segmentation model */
    SegmentationModel model(alpha_prefix, alpha_stem, alpha_suffix,
           word_vocabulary, substring_vocabulary.Size(), tries);

    std::random_device rd;
    std::mt19937 engine(rd());

    std::vector<Segmentation> segs;

    /* Obtain initial random segmentations */
    unsigned wid = 0;
    for(auto& sentence: corpus) {
        for(auto& word: sentence) {
            const Segmentation seg = model.Increment(word, engine, true);
            segs.push_back(seg);
            wid++;
        }
    }

    std::cerr << "Initialization done \n"
        << "Running parallel Gibbs sampler with " << NTHREADS << " threads\n";

    /* Run Gibbs sampler */
    for(unsigned it = 0; it < n_iterations; it++) {
        wid = 0;
        ThreadPool pool(NTHREADS);
        for(auto& sentence: corpus) {
            for(auto word: sentence) {
                pool.enqueue([&model, &engine, &segs, wid, word] {
                model.Decrement(word, segs[wid]);
                segs[wid] = model.Increment(word, engine, false);
                });
                wid++;
            }
        }
        pool.join();

        if(it % 10 == 0) {
            std::cerr << "Iteration " << (it+1) << "/" << n_iterations << "\n";
            std::cerr << model << "\n";
            double ll = model.LogLikelihood();
            double ppl = exp(-ll/corpus.Tokens());
            std::cerr << "LL=" << ll << " ppl=" << ppl << "\n";
        }
    }

    /* Print final segmentations decoded with Viterbi algorithm */
    for(unsigned w = 0; w < word_vocabulary.Size(); w++) {
        const std::string& word = word_vocabulary.Convert(w);
        const Segmentation seg = model.Decode(w);
        std::cout << word << "\t" << FormatSegmentation(seg, substring_vocabulary) << "\n";
    }
}
