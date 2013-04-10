#include <map>
#include <iostream>

struct Trie {
    std::map<char, Trie> nodes;
    int label;

    Trie() : nodes(), label(-1) {}

    void Insert(const std::string& str, int label) {
        if(str.size() == 0)
            this->label = label;
        else
            nodes[str[0]].Insert(str.substr(1), label);
    }

    void Print(std::ostream& out, const std::string& prefix="") {
        if(label != -1)
            out << prefix << " -> " << label << "\n";
        for(auto& child: nodes) {
            child.second.Print(out, prefix+child.first);
        }
    }
};
