struct Segmentation {
    vector<unsigned> prefixes, suffixes;
    unsigned stem;
};


template <typename Arc>
const Segmentation ReadSegmentation(const fst::ExpandedFst<Arc>& path, const Trie& trie) {
    vector<unsigned> prefixes, suffixes;
    unsigned stem = -1;
    unsigned part = 0;
    const Trie* node = &trie;
    for(fst::StateIterator<fst::ExpandedFst<Arc>> siter(path); !siter.Done(); siter.Next()) {
        typename fst::ExpandedFst<Arc>::StateId state_id = siter.Value();
        for(fst::ArcIterator<fst::ExpandedFst<Arc>> aiter(path, state_id);
                !aiter.Done(); aiter.Next()) {
            const Arc &arc = aiter.Value();
            if(arc.olabel == '^') {
                (part == 0 ? prefixes : suffixes).push_back(node->label);
                //std::cerr << "^" << node->label << "^";
                node = &trie;
            }
            else if(arc.olabel == '<') {
                //std::cerr << "<";
                part++;
            }
            else if(arc.olabel == '>') {
                stem = node->label;
                //std::cerr << node->label << ">";
                node = &trie;
                part++;
            }
            else node = &node->nodes.at(arc.olabel);
        }
    }
    //std::cerr << "\n";
    assert(stem != -1);
    return Segmentation {prefixes, suffixes, stem};
}

class SegmentationModel {
    public:
    SegmentationModel(float alpha_prefix, float alpha_stem, float alpha_suffix,
            const Vocabulary& word_vocabulary, unsigned n_substrings,
            const std::vector<Trie>& tries) :
        word_vocabulary(word_vocabulary),
        tries(tries),
        prefix_model(n_substrings, alpha_prefix),
        stem_model(n_substrings, alpha_stem),
        suffix_model(n_substrings, alpha_suffix),
        prefix_length_model(1, 1),
        suffix_length_model(1, 1) {}

    const Segmentation Increment(unsigned w,
            std::mt19937& engine, bool initialize=false) {
        const std::string& word = word_vocabulary.Convert(w);
        fst::LogVectorFst log_lattice = MakeLattice<fst::LogArc>(tries[w], word);
        fst::LogVectorFst sampled;
        int seed = prob::randint(engine, -INT_MAX, INT_MAX);
        if(initialize) {
            fst::UniformArcSelector<fst::LogArc> selector(seed);
            fst::RandGenOptions< fst::UniformArcSelector<fst::LogArc> > options(selector);
            fst::RandGen(log_lattice, &sampled);
        }
        else {
            std::vector<fst::LogWeight> beta;
            fst::ShortestDistance<fst::LogArc>(log_lattice, &beta, true);
            fst::Reweight<fst::LogArc>(&log_lattice, beta, fst::REWEIGHT_TO_INITIAL);
            fst::LogProbArcSelector<fst::LogArc> selector(seed);
            fst::RandGenOptions< fst::LogProbArcSelector<fst::LogArc> > options(selector);
            fst::RandGen(log_lattice, &sampled, options);
        }

        const Segmentation seg = ReadSegmentation(sampled, tries[w]);
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
        const fst::StdVectorFst lattice = MakeLattice<fst::StdArc>(tries[w], word);
        fst::StdVectorFst best;
        fst::ShortestPath(lattice, &best);
        fst::TopSort(&best);
        const Segmentation seg = ReadSegmentation(best, tries[w]);
        return seg;
    }

    double LogLikelihood() const {
        return prefix_model.LogLikelihood()
            + stem_model.LogLikelihood()
            + suffix_model.LogLikelihood();
        // TODO + length likelhood
    }

    private:
    template <typename Arc>
    fst::VectorFst<Arc> MakeLattice(const Trie& trie, const std::string& word) const {
        const fst::VectorFst<Arc> grammar = BuildGrammar<Arc>(trie, 
                prefix_model, stem_model, suffix_model,
                prefix_length_model, suffix_length_model);
        //grammar.Write("grammar.fst");

        const fst::VectorFst<Arc> word_fst = LinearChain<Arc>(word);
        //word_fst.Write("chain.fst");

        fst::VectorFst<Arc> lattice;
        fst::Compose(word_fst, grammar, &lattice);
        fst::RmEpsilon(&lattice);
        //lattice.Write("lattice.fst");

        return lattice;
    }

    const Vocabulary& word_vocabulary;
    const std::vector<Trie>& tries;
    DirichletMultinomial prefix_model, stem_model, suffix_model;
    BetaGeometric prefix_length_model, suffix_length_model;

    friend std::ostream& operator<<(std::ostream&, const SegmentationModel&);
};


std::ostream& operator<<(std::ostream& os, const SegmentationModel& m) {
    return os << "SegmentationModel(prefix ~ " << m.prefix_model
        << ", |prefix| ~ " << m.prefix_length_model
        << ", stem ~ " << m.stem_model
        << ", suffix ~ " << m.suffix_model
        << ", |suffix| ~ " << m.suffix_length_model << ")";
}
