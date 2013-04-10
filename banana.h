namespace fst {
typedef VectorFst<LogArc> LogVectorFst;
}

void BuildBanana(const Trie& trie, int start, int end,
        fst::LogVectorFst &grammar, const DirichletMultinomial& model) {
    for(auto& child: trie.nodes) {
        char c = child.first;
        int label = child.second.label;
        float w = -log(model.Prob(label));
        // ilabel, olabel, weight, dest state ID.
        grammar.AddArc(start, fst::LogArc((unsigned char) c, (unsigned char) c, w, end));
        if(child.second.nodes.size() > 0) {
            const int k = grammar.AddState();
            grammar.AddArc(start, fst::LogArc((unsigned char) c, (unsigned char) c, 0, k));
            BuildBanana(child.second, k, end, grammar, model);
        }
    }
}

const int mb = (unsigned char)'^';
const int ss = (unsigned char)'<';
const int se = (unsigned char)'>';

const fst::LogVectorFst BuildGrammar(const Trie& trie,
        const DirichletMultinomial& prefix_model,
        const DirichletMultinomial& stem_model,
        const DirichletMultinomial& suffix_model,
        const BetaGeometric& prefix_length_model,
        const BetaGeometric& suffix_length_model) {
    fst::LogVectorFst grammar;

    // Prefix
    const float prefix_loop = -log(1 - prefix_length_model.Stop());
    const float prefix_stop = -log(prefix_length_model.Stop());
    const int prefix_start = grammar.AddState();
    grammar.SetStart(prefix_start);
    const int prefix_end = grammar.AddState();
    const int prefix_m = grammar.AddState();
    grammar.AddArc(prefix_end, fst::LogArc(0, mb, prefix_loop, prefix_m)); // morpheme penalty
    grammar.AddArc(prefix_m, fst::LogArc(0, 0, 0, prefix_start)); // closure
    grammar.AddArc(prefix_start, fst::LogArc(0, 0, 0, prefix_m)); // closure
    BuildBanana(trie, prefix_start, prefix_end, grammar, prefix_model);

    // Stem
    const int stem_start = grammar.AddState();
    grammar.AddArc(prefix_m, fst::LogArc(0, ss, prefix_stop, stem_start)); // start stem
    const int stem_end = grammar.AddState();
    BuildBanana(trie, stem_start, stem_end, grammar, stem_model);

    // Suffix
    const float suffix_loop = -log(1 - suffix_length_model.Stop());
    const float suffix_stop = -log(suffix_length_model.Stop());
    const int suffix_start = grammar.AddState();
    grammar.AddArc(stem_end, fst::LogArc(0, se, 0, suffix_start)); // end stem
    const int suffix_end = grammar.AddState();
    const int suffix_m = grammar.AddState();
    grammar.AddArc(suffix_end, fst::LogArc(0, mb, suffix_loop, suffix_m)); // morpheme penalty
    grammar.AddArc(suffix_m, fst::LogArc(0, 0, 0, suffix_start)); // closure
    grammar.AddArc(suffix_start, fst::LogArc(0, 0, 0, suffix_m)); // closure
    BuildBanana(trie, suffix_start, suffix_end, grammar, suffix_model);
    const int final = grammar.AddState();
    grammar.AddArc(suffix_m, fst::LogArc(0, 0, suffix_stop, final)); // end word
    grammar.SetFinal(final, 0);

    return grammar;
}

const fst::LogVectorFst LinearChain(const std::string& word) {
    fst::LogVectorFst chain;
    for(unsigned i = 0; i < word.size(); i++) {
        chain.AddState();
        chain.AddArc(i, fst::LogArc((unsigned char) word[i], (unsigned char) word[i], 0, i+1));
        //std::cerr << (unsigned char) word[i] << ", ";
    }
    //std::cerr << "\n";
    chain.SetStart(0);
    chain.SetFinal(chain.AddState(), 0);
    return chain;
}
