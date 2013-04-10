struct Segmentation {
    vector<unsigned> prefixes, suffixes;
    unsigned stem;
};


template <typename Arc>
const Segmentation ReadSegmentation(const fst::ExpandedFst<Arc>& path,
        const Vocabulary& substring_vocabulary) {
    assert(path.NumStates() > 0);
    vector<unsigned> prefixes, suffixes;
    unsigned stem = substring_vocabulary.Size();
    unsigned part = 0;
    std::string morpheme; // TODO use trie instead
    for(fst::StateIterator<fst::ExpandedFst<Arc>> siter(path); !siter.Done(); siter.Next()) {
        typename fst::ExpandedFst<Arc>::StateId state_id = siter.Value();
        for(fst::ArcIterator<fst::ExpandedFst<Arc>> aiter(path, state_id);
                !aiter.Done(); aiter.Next()) {
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
    //std::cerr << " / " << stem << " < " << substring_vocabulary.Size() << "\n";
    assert(stem < substring_vocabulary.Size());
    return Segmentation {prefixes, suffixes, stem};
}

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

    const Segmentation Increment(unsigned w,
            std::mt19937& engine, bool initialize=false) {
        const std::string& word = word_vocabulary.Convert(w);
        fst::LogVectorFst lattice = MakeLattice(tries[w], word);
        fst::LogVectorFst sampled;
        if(initialize) {
            fst::UniformArcSelector<fst::LogArc> selector; // FIXME seed
            fst::RandGenOptions< fst::UniformArcSelector<fst::LogArc> > options(selector);
            fst::RandGen(lattice, &sampled);
        }
        else {
            std::vector<fst::LogWeight> beta;
            fst::ShortestDistance<fst::LogArc>(lattice, &beta, true);
            fst::Reweight<fst::LogArc>(&lattice, beta, fst::REWEIGHT_TO_INITIAL);
            fst::LogProbArcSelector<fst::LogArc> selector; // FIXME seed
            fst::RandGenOptions< fst::LogProbArcSelector<fst::LogArc> > options(selector);
            fst::RandGen(lattice, &sampled, options);
        }

        const Segmentation seg = ReadSegmentation(sampled, substring_vocabulary);
        for(unsigned p: seg.prefixes)
            prefix_model.Increment(p);
        prefix_length_model.Increment(seg.prefixes.size());
        stem_model.Increment(seg.stem);
        for(unsigned s: seg.suffixes)
            suffix_model.Increment(s);
        suffix_length_model.Increment(seg.suffixes.size());
        return seg;
    }

    void Decrement(unsigned w, const Segmentation& seg) {
        const std::string& word = word_vocabulary.Convert(w);
        for(unsigned p: seg.prefixes)
            prefix_model.Decrement(p);
        prefix_length_model.Decrement(seg.prefixes.size());
        stem_model.Decrement(seg.stem);
        for(unsigned s: seg.suffixes)
            suffix_model.Decrement(s);
        suffix_length_model.Decrement(seg.suffixes.size());
    }

    const Segmentation Decode(unsigned w) const {
        const std::string& word = word_vocabulary.Convert(w);
        const fst::LogVectorFst log_lattice = MakeLattice(tries[w], word);
        fst::StdVectorFst lattice;
        fst::WeightConvertMapper<fst::LogArc, fst::StdArc> mapper;
        fst::ArcMap(log_lattice, &lattice, mapper);
        fst::StdVectorFst best;
        fst::ShortestPath(lattice, &best);
        fst::TopSort(&best);
        const Segmentation seg = ReadSegmentation(best, substring_vocabulary);
        return seg;
    }

    double LogLikelihood() const {
        return prefix_model.LogLikelihood()
            + stem_model.LogLikelihood()
            + suffix_model.LogLikelihood();
        // TODO + length likelhood
    }

    private:
    // FIXME templatize by arc type (as well as BuildGrammar, LC...)
    fst::LogVectorFst MakeLattice(const Trie& trie, const std::string& word) const {
        //std::cerr << word << "\n";
        const fst::LogVectorFst grammar = BuildGrammar(trie, 
                prefix_model, stem_model, suffix_model,
                prefix_length_model, suffix_length_model);
        //std::cerr << grammar.NumStates() << "\n";

        const fst::LogVectorFst word_fst = LinearChain(word);
        //std::cerr << word_fst.NumStates() << "\n";

        fst::LogVectorFst lattice;
        fst::Compose(word_fst, grammar, &lattice);
        fst::RmEpsilon(&lattice);
        fst::Project(&lattice, fst::PROJECT_OUTPUT);
        //std::cerr << lattice.NumStates() << "\n";

        if(lattice.NumStates() == 0) {
            grammar.Write("grammar.fst");
            word_fst.Write("chain.fst");
        }

        return lattice;
    }

    const Vocabulary& word_vocabulary;
    const Vocabulary& substring_vocabulary;
    const std::vector<Trie>& tries;
    DirichletMultinomial prefix_model, stem_model, suffix_model;
    BetaGeometric prefix_length_model, suffix_length_model;
};
