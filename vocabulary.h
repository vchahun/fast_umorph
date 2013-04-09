#include <map>
#include <vector>
#include <string>

class Vocabulary {
    std::map<std::string, unsigned> word2id;
    std::vector<std::string> id2word;

    public:
    Vocabulary() : word2id(), id2word() {}

    unsigned Encode(const std::string& word) {
        std::map<std::string, unsigned>::const_iterator it = word2id.find(word);
        if(it != word2id.end())
            return it->second;
        id2word.push_back(word);
        word2id[word] = id2word.size()-1;
        return id2word.size()-1;
    }

    unsigned Convert(const std::string& word) const {
        return word2id.at(word);
    }

    const std::string& Convert(unsigned k) const {
        assert(k < id2word.size());
        return id2word[k];
    }

    size_t Size() const {
        return id2word.size();
    }

};
