function list = discreteTimeSeries2Words(words, nAlphabet)
% list = discreteTimeSeries2Words(words, nAlphabet)
% Finds the unique words and replaces each word with a corresponding integer
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
%   counts: (wordLength x 1) integer representation of the same word sequence
%
% See also: words2multiplicities, testFastWords, multiplicitiesFromCounts, fastWord2Counts
%
% To compile the C-MEX file:
%   mex discreteTimeSeries2Words.c fastWordsTree.c
