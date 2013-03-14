%% Example code from the README file.

%% Convert to multiplicities from samples, using multiplicitiesFromSamples
[mm, icts] = multiplicitiesFromSamples([1 2 2 3.5 3.5 4 4 4 4]);
assert(all(mm==[1 2 1]'))
assert(all(icts==[1 2 4]'))

%% Convert to multiplicities from samples, using multiplicitiesFromSamples (uint16) 
[mm, icts] = multiplicitiesFromSamples(uint16([1 2 2 3 3 4 4 4 4]));
assert(all(mm==[1 2 1]'))
assert(all(icts==[1 2 4]'))

%% Convert to multiplicities from histogram representation
[mm, icts] = multiplicitiesFromCounts([1 2 2 4]);
assert(all(mm==[1 2 1]'))
assert(all(icts==[1 2 4]'))

%% Convert to histogram representation from multiplicities.
mm = [1 2 1]; icts = [1 2 4];
hg = multiplicitiesToCounts(mm, icts);
assert(all(hg == [4 2 2 1]'))
fprintf('\nExample code successful!\n')

