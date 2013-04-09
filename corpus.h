#include <string>
#include <vector>
#include <iostream>
#include <sstream>

class Corpus {
    std::vector< std::vector<unsigned> > segments;
    typedef std::vector< std::vector<unsigned> >::const_iterator const_iterator;
    
    public:
    Corpus(std::istream& input_stream, Vocabulary& vocabulary): segments() {
        std::string line, word;
        while(getline(input_stream, line)) {
            std::vector<unsigned> segment;
            std::stringstream sstr(line);
            while(sstr >> word) {
                segment.push_back(vocabulary.Encode(word));
            }
            segments.push_back(segment);
        }
    }

    size_t Size() const {
        return segments.size();
    }

    size_t Tokens() const {
        size_t N = 0;
        for(auto& doc: segments) {
            N += doc.size();
        }
        return N;
    }

    const_iterator begin() const {
        return segments.begin(); 
    }

    const_iterator end() const {
        return segments.end(); 
    }
};

