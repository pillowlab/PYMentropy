wordLength = 8;
nWords = 100;
words = uint16(zeros(wordLength, nWords));

words(3, 10) = 1;
words(4, 10) = 1;

words(4, 11) = 1;
words(5, 11) = 1;

words(4, 20) = 1;
words(4, 30) = 1;
words(4, 40) = 1;

h = fastWords2Counts(words, 2)
assert(all(h == [1;1;3;95]));
[mm, icts] = multiplicitiesFromCounts(h)
[mm, icts] = words2multiplicities(double(words)')

x = discreteTimeSeries2Words(words, 2);
assert(length(unique(x)) == 4);
