%% DON'T CHANGE THIS UNLESS YOU WANT TO OVERWRITE THE CURRENT TEST VALUES!
GENERATE_TEST = 0; % set = 1 to generate a new set of test values

if(GENERATE_TEST)
    Ntest = 25;
    seedval = 12123;
    H_TEST = zeros(Ntest,2);
else
    tstr = load('test_Hvals.mat');
    Ntest = tstr.Ntest;
    seedval = tstr.seedval;
    H_TEST = tstr.H_TEST;
end

assertRange = @(x,y) assert(abs(x-y) < 2e-11, 'err: [%g]', x-y);

rand('seed', seedval)
randn('seed', seedval)
% s = RandStream('mcg16807','Seed', seedval);
% RandStream.setDefaultStream(s);

param = struct('nagrid', 25, 'ndgrid', 25);

for idx = 1:Ntest
    fprintf('Running test %d\r', idx)
    K = round(50*rand+5);
    N = round(50*rand+5);
    nn = zeros(K,1);
    while( sum(nn) < 2 || length(icts) == 1)
        nn = mnrnd(N,dirichletrnd(2*rand, K));
        [mm icts] = multiplicitiesFromCounts(nn);
    end
    [Hbls, Hvar] = computeH_PYM(mm, icts, param);
    if(GENERATE_TEST)
        H_TEST(idx,1) = Hbls; H_TEST(idx,2) = Hvar;
    else
        assertRange(Hbls, H_TEST(idx,1));
        if(Hvar==inf)
            assert(Hvar == H_TEST(idx,2))
        else
            assertRange(Hvar, H_TEST(idx,2));
        end
        fprintf('\t.... passed!\n')
    end
end

if(GENERATE_TEST)
    save('test_Hvals.mat', 'H_TEST', 'seedval', 'Ntest')
end

%% Finally, some hand-built tests in case all this random nonsense falls apart.

% [mm1, icts1] = multiplicitiesFromCounts([2 1 1 1]);
% [mm2, icts2] = multiplicitiesFromCounts([2 2 1 1]);
% [mm3, icts3] = multiplicitiesFromCounts([5 3 2 1]);
% 
% [Hbls, Hvar] = computeH_PYM(mm1, icts1);
% assertRange(Hbls, 3.965631191137517); assert(Hvar == inf);
% 
% [Hbls, Hvar] = computeH_PYM(mm2, icts2);
% assertRange(Hbls, 3.067407634404604);
% assertRange(Hvar, 0.522879923326087);
% 
% [Hbls, Hvar] = computeH_PYM(mm3, icts3);
% assertRange(Hbls, 1.894590969625430);
% assertRange(Hvar, 0.195283756259754);
% 
% fprintf('\tPassed hand-coded tests!\n')


% K = round(50*rand);
% N = round(50*rand);
% nn1 = rand(K,1); nn1 = nn1/sum(nn1); nn1 = mnrnd(N,nn1);
% [mm1 icts1] = multiplicitiesFromCounts(nn1);


