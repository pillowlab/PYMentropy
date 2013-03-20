function counts = fastWords2Counts(words, nAlphabet)
% counts = fastWords2Counts(words, nAlphabet)
% Creates a histogram of word frequencies from a matrix representation of words
%
% Uses a nAlphabet-ary tree that only instantiates for the occurring words.
%
% Input
%   words: (wordLength x nWords; uint32) 
%	each word is a column vector of length wordLength.
%	each word is composed of alphabets of non-negative integers < nAlphabet
%   nAlphabet: each alphabet is an integer in [0, 1, ..., nAlphabet-1].
%	use the smallest nAlphabet for memory and time efficiency
%
% Output
%   counts: (M x 1) histogram of word occurrence frequency
%
% See also: words2multiplicities, testFastWords, multiplicitiesFromCounts, discreteTimeSeries2Words
%
% To compile the C-MEX file:
%   mex fastWords2Counts.c fastWordsTree.c
