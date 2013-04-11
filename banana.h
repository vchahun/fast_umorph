namespace fst {
typedef VectorFst<LogArc> LogVectorFst;
}

// Add arcs corresponding to 'trie' between 'start' and 'end' nodes
// For each trie leaf, add the corresponding weight in 'model'
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
    const int prefix1 = grammar.AddState(); // start -> 1 (closure)
    grammar.AddArc(prefix_start, Arc(0, 0, 0, prefix1));
    const int prefix2 = grammar.AddState(); // 1 -> substrings -> 2
    BuildBanana(trie, prefix1, prefix2, grammar, prefix_model);
    const int prefix3 = grammar.AddState(); // 2 -> 3 / p; 3 -> 1 (closure)
    grammar.AddArc(prefix2, Arc(0, mb, prefix_loop, prefix3)); // morpheme penalty
    grammar.AddArc(prefix3, Arc(0, 0, 0, prefix1)); // closure
    const int prefix_end = grammar.AddState(); // start, 3 -> end (closure)
    grammar.AddArc(prefix_start, Arc(0, 0, 0, prefix_end)); // closure
    grammar.AddArc(prefix3, Arc(0, 0, 0, prefix_end)); // closure

    // Stem
    const int stem_start = grammar.AddState();
    grammar.AddArc(prefix_end, Arc(0, ss, prefix_stop, stem_start)); // prefix -> suffix
    const int stem_end = grammar.AddState();
    BuildBanana(trie, stem_start, stem_end, grammar, stem_model);

    // Suffix
    const float suffix_loop = -log(1 - suffix_length_model.Stop());
    const float suffix_stop = -log(suffix_length_model.Stop());
    const int suffix_start = grammar.AddState();
    grammar.AddArc(stem_end, Arc(0, se, 0, suffix_start)); // stem -> suffix
    const int suffix1 = grammar.AddState(); // start -> 1 (closure)
    grammar.AddArc(suffix_start, Arc(0, 0, 0, suffix1));
    const int suffix2 = grammar.AddState(); // 1 -> substrings -> 2
    BuildBanana(trie, suffix1, suffix2, grammar, suffix_model);
    const int suffix3 = grammar.AddState(); // 2 -> 3 / p; 3 -> 1 (closure)
    grammar.AddArc(suffix2, Arc(0, mb, suffix_loop, suffix3)); // morpheme penalty
    grammar.AddArc(suffix3, Arc(0, 0, 0, suffix1)); // closure
    const int suffix_end = grammar.AddState(); // start, 3 -> end (closure)
    grammar.AddArc(suffix_start, Arc(0, 0, 0, suffix_end)); // closure
    grammar.AddArc(suffix3, Arc(0, 0, 0, suffix_end)); // closure
    grammar.SetFinal(suffix_end, suffix_stop);

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
