# A simple Bayesian word segmentation model

Code for unsupervised morphology/segmentation induction.

## Model

We assume a grammar of the form M\*MM\* where M is the set of all possible substrings. Each word should have possibly one or more prefixes, a stem and possibly one or more suffixes.

The following generative story is used:

    theta_prefix ~ Dir(alpha_prefix)
    theta_stem ~ Dir(alpha_stem)
    theta_suffix ~ Dir(alpha_suffix)

    lambda_prefix ~ Uniform([0, 1])
    lambda_suffix ~ Uniform([0, 1])

    for each word:
        sample lenghts |prefixes| ~ Geometric(lambda_prefix),
            |suffixes| ~ Geometric(lambda_suffix)
        sample prefixes, a stem and suffixes from the theta_. distributions
        concatenate prefixes, the stem and suffixes to form the word

The likelihood decomposes into factors representable by a WFSA.
We use the OpenFst library to implement Gibbs sampling with this model.

## Compiling

The code requires a C++11 compatible compiler (the Makefile assumes gcc 4.7) and the [OpenFst library](http://www.openfst.org/) has to be installed.

Execute `make` to compile the `segment` binary.

## Running

Create a tokenized corpus or a vocabulary file and pipe it into the program:

    cat words.txt | ./segment 1000 1e-5 1e-6 1e-5 > words.segs.txt

The corpus has to be encoded in a **8-bit format** (e.g. ISO-8859-1 for French, ISO-8859-5 for Russian, ISO-8859-6 for Arabic...), due to my unforgivable laziness to properly treat unicode strings. You might want to run a pipeline similar to:

    cat french-words.txt | iconv -f utf8 -t latin1 |\
        ./segment 1000 1e-5 1e-6 1e-5 |\
        iconv -f latin1 -t utf8 > french-words.segs.txt

Also, the characters `^<>` are currently reserved as special morpheme boundary markers but this can easily be changed in the code.

## Parameters

The parameters `alpha_prefix`, `alpha_stem` and `alpha_suffix` have to be adjusted to force segmentation. In practice, we have found that `alpha_prefix, alpha_suffix << alpha_stem << 1` is necessary to obtain useful segmentations for a variety of languages.

## License

Copyright (c) 2013, [Victor Chahuneau](http://victor.chahuneau.fr/)

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
