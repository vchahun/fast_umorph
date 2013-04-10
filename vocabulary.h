#include <unordered_map>
#include <vector>
#include <string>
#include <cassert>

class Vocabulary {
    std::unordered_map<std::string, unsigned> word2id;
    std::vector<std::string> id2word;
    typedef std::vector<std::string>::const_iterator const_iterator;

    public:
    Vocabulary() : word2id(), id2word() {}

    unsigned Encode(const std::string& word) {
        std::unordered_map<std::string, unsigned>::const_iterator it = word2id.find(word);
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

    const_iterator begin() const {
        return id2word.begin(); 
    }

    const_iterator end() const {
        return id2word.end(); 
    }

};
