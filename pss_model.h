template <typename Arc>
struct Segmentation {
    vector<unsigned> prefixes;
    vector<unsigned> suffixes;
    unsigned stem;

    typedef typename Arc::StateId StateId;
    typedef fst::Fst<Arc> ArcFst;
    typedef fst::StateIterator<ArcFst> state_iterator;
    typedef fst::ArcIterator<ArcFst> arc_iterator;

    Segmentation(const ArcFst& path, const Vocabulary& substring_vocabulary) {
        unsigned part = 0;
        std::string morpheme; // TODO use trie instead
        for(state_iterator siter(path); !siter.Done(); siter.Next()) {
            StateId state_id = siter.Value();
            for(arc_iterator aiter(path, state_id); !aiter.Done(); aiter.Next()) {
                const Arc &arc = aiter.Value();
                if(arc.olabel == '^') {
                    //std::cerr << "^" << morpheme << "^";
                    unsigned m = substring_vocabulary.Convert(morpheme);
                    (part == 0 ? prefixes : suffixes).push_back(m);
                    morpheme = "";
                }
                else if(arc.olabel == '<') {
                    part++;
                    //std::cerr << "<";
                }
                else if(arc.olabel == '>') {
                    stem = substring_vocabulary.Convert(morpheme);
                    //std::cerr << morpheme << ">";
                    morpheme = "";
                    part++;
                }
                else morpheme += arc.olabel;
            }
        }
        //std::cerr << "\n";
    }
};

class SegmentationModel {
    public:
    SegmentationModel(float alpha_prefix, float alpha_stem, float alpha_suffix,
            const Vocabulary& word_vocabulary, const Vocabulary& substring_vocabulary,
            const std::vector<Trie>& tries) :
        word_vocabulary(word_vocabulary),
        substring_vocabulary(substring_vocabulary),
        tries(tries),
        prefix_model(substring_vocabulary.Size(), alpha_prefix),
        stem_model(substring_vocabulary.Size(), alpha_stem),
        suffix_model(substring_vocabulary.Size(), alpha_suffix),
        prefix_length_model(1, 1),
        suffix_length_model(1, 1) {}

    const Segmentation<fst::LogArc> Increment(unsigned w,
            std::mt19937& engine, bool initialize=false) {
        const std::string& word = word_vocabulary.Convert(w);
        const fst::StdVectorFst lat = MakeLattice(tries[w], word);
        fst::VectorFst<fst::LogArc> lattice;
        fst::WeightConvertMapper<fst::StdArc, fst::LogArc> mapper;
        fst::ArcMap(lat, &lattice, mapper);
        fst::VectorFst<fst::LogArc> sampled;
        if(initialize) {
            fst::RandGen(lattice, &sampled);
        }
        else {
            // FIXME original FST should have log weights? / convert
            std::vector<fst::LogWeight> beta;
            fst::ShortestDistance<fst::LogArc>(lattice, &beta, true);
            fst::Reweight(&lattice, beta, fst::REWEIGHT_TO_INITIAL);
            fst::LogProbArcSelector<fst::LogArc> selector;
            fst::RandGenOptions< fst::LogProbArcSelector<fst::LogArc> > options(selector);
            fst::RandGen(lattice, &sampled, options);
        }

        Segmentation<fst::LogArc> seg(sampled, substring_vocabulary);
        for(unsigned p: seg.prefixes)
            prefix_model.Increment(p);
        prefix_length_model.Increment(seg.prefixes.size());
        stem_model.Increment(seg.stem);
        for(unsigned s: seg.suffixes)
            suffix_model.Increment(s);
        suffix_length_model.Increment(seg.suffixes.size());
        return seg;
    }

    void Decrement(unsigned w, const Segmentation<fst::LogArc>& seg) {
        const std::string& word = word_vocabulary.Convert(w);
        for(unsigned p: seg.prefixes)
            prefix_model.Decrement(p);
        prefix_length_model.Decrement(seg.prefixes.size());
        stem_model.Decrement(seg.stem);
        for(unsigned s: seg.suffixes)
            suffix_model.Decrement(s);
        suffix_length_model.Decrement(seg.suffixes.size());
    }

    const Segmentation<fst::StdArc> Decode(unsigned w) const {
        const std::string& word = word_vocabulary.Convert(w);
        fst::StdVectorFst best;
        fst::ShortestPath(MakeLattice(tries[w], word), &best);
        fst::TopSort(&best);
        Segmentation<fst::StdArc> seg(best, substring_vocabulary);
        return seg;
    }

    double LogLikelihood() const {
        return prefix_model.LogLikelihood()
            + stem_model.LogLikelihood()
            + suffix_model.LogLikelihood();
        // TODO + length likelhood
    }

    private:
    const fst::StdVectorFst MakeLattice(const Trie& trie, const std::string& word) const {
        const fst::StdVectorFst grammar = BuildGrammar(trie, 
                prefix_model, stem_model, suffix_model,
                prefix_length_model, suffix_length_model);

        const fst::StdVectorFst word_fst = LinearChain(word);

        fst::StdVectorFst lattice;
        fst::Compose(word_fst, grammar, &lattice);
        fst::RmEpsilon(&lattice);
        fst::Project(&lattice, fst::PROJECT_OUTPUT);
        return lattice;
    }

    const Vocabulary& word_vocabulary;
    const Vocabulary& substring_vocabulary;
    const std::vector<Trie>& tries;
    DirichletMultinomial prefix_model, stem_model, suffix_model;
    BetaGeometric prefix_length_model, suffix_length_model;
};
