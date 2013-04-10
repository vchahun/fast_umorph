namespace fst {
typedef VectorFst<LogArc> LogVectorFst;
}

template <typename Arc>
void BuildBanana(const Trie& trie, int start, int end,
        fst::VectorFst<Arc> &grammar, const DirichletMultinomial& model) {
    for(auto& child: trie.nodes) {
        char c = child.first;
        int label = child.second.label;
        float w = -log(model.Prob(label));
        grammar.AddArc(start, Arc((unsigned char) c, (unsigned char) c, w, end));
        if(child.second.nodes.size() > 0) {
            const int k = grammar.AddState();
            grammar.AddArc(start, Arc((unsigned char) c, (unsigned char) c, 0, k));
            BuildBanana(child.second, k, end, grammar, model);
        }
    }
}

const int mb = (unsigned char)'^';
const int ss = (unsigned char)'<';
const int se = (unsigned char)'>';

template <typename Arc>
const fst::VectorFst<Arc> BuildGrammar(const Trie& trie,
        const DirichletMultinomial& prefix_model,
        const DirichletMultinomial& stem_model,
        const DirichletMultinomial& suffix_model,
        const BetaGeometric& prefix_length_model,
        const BetaGeometric& suffix_length_model) {
    fst::VectorFst<Arc> grammar;

    // Prefix
    const float prefix_loop = -log(1 - prefix_length_model.Stop());
    const float prefix_stop = -log(prefix_length_model.Stop());
    const int prefix_start = grammar.AddState();
    grammar.SetStart(prefix_start);
    const int prefix_end = grammar.AddState();
    const int prefix_m = grammar.AddState();
    grammar.AddArc(prefix_end, Arc(0, mb, prefix_loop, prefix_m)); // morpheme penalty
    grammar.AddArc(prefix_m, Arc(0, 0, 0, prefix_start)); // closure
    grammar.AddArc(prefix_start, Arc(0, 0, 0, prefix_m)); // closure
    BuildBanana(trie, prefix_start, prefix_end, grammar, prefix_model);

    // Stem
    const int stem_start = grammar.AddState();
    grammar.AddArc(prefix_m, Arc(0, ss, prefix_stop, stem_start)); // start stem
    const int stem_end = grammar.AddState();
    BuildBanana(trie, stem_start, stem_end, grammar, stem_model);

    // Suffix
    const float suffix_loop = -log(1 - suffix_length_model.Stop());
    const float suffix_stop = -log(suffix_length_model.Stop());
    const int suffix_start = grammar.AddState();
    grammar.AddArc(stem_end, Arc(0, se, 0, suffix_start)); // end stem
    const int suffix_end = grammar.AddState();
    const int suffix_m = grammar.AddState();
    grammar.AddArc(suffix_end, Arc(0, mb, suffix_loop, suffix_m)); // morpheme penalty
    grammar.AddArc(suffix_m, Arc(0, 0, 0, suffix_start)); // closure
    grammar.AddArc(suffix_start, Arc(0, 0, 0, suffix_m)); // closure
    BuildBanana(trie, suffix_start, suffix_end, grammar, suffix_model);
    const int final = grammar.AddState();
    grammar.AddArc(suffix_m, Arc(0, 0, suffix_stop, final)); // end word
    grammar.SetFinal(final, 0);

    return grammar;
}

template <typename Arc>
const fst::VectorFst<Arc> LinearChain(const std::string& word) {
    fst::VectorFst<Arc> chain;
    for(unsigned i = 0; i < word.size(); i++) {
        chain.AddState();
        chain.AddArc(i, Arc((unsigned char) word[i], (unsigned char) word[i], 0, i+1));
    }
    chain.SetStart(0);
    chain.SetFinal(chain.AddState(), 0);
    return chain;
}
