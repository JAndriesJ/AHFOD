%% For each night pair, form the random dist
function [RandomDistofSP, ActualDistofSP] = getTestReTest(HFOrateMat, nbPermutations)
    if nargin < 2
    nbPermutations = 100000;
    end
    
    nbRecordings     = size(HFOrateMat,1);
    rateMatrix       = HFOrateMat'; % channel X Interval
    [RandomDistofSP, ActualDistofSP] = getPairswiseRandomSP(rateMatrix, nbRecordings, nbPermutations);
end

%% Getting the random distribution:
function [SortedSPAllNightPairs, pairswiseSP] = getPairswiseRandomSP(rateMatrix, nbRecordings, M )
    AllNightCombinations = nchoosek(1:nbRecordings,2);
    nbNightCombos = size(AllNightCombinations,1);
    SortedSPAllNightPairs = nan(M/2, nbNightCombos);
    pairswiseSP = nan(1,nbNightCombos);
    for iComb = 1:nbNightCombos
        InterVal1 = AllNightCombinations(iComb,1);
        InterVal2 = AllNightCombinations(iComb,2);

        RatePair = rateMatrix(:, [InterVal1, InterVal2]);
        
        SortedSP = getRandomDistribution(RatePair, M);
        SortedSPAllNightPairs(:,iComb) = SortedSP(:);
        pairswiseSP(iComb) =  (RatePair(:,1)/norm(RatePair(:,1)))'*(RatePair(:,2)/norm(RatePair(:,2)));
    end
end

function [SortedScalarProduct ] = getRandomDistribution(RatePair, M )
    N = size(RatePair,1);
    PermutationBank = Get_M_Combinations_to_N(M,N);
    % M = min(size(PermutationBank,1),M/2);
    CombinationsForNight1 = PermutationBank(1:M/2,:);
    CombinationsForNight2 = PermutationBank(M/2+1:M,:);
    RatesForNights1 = RatePair(:,1);
    RatesForNights2 = RatePair(:,2);

    ScalarProduct = zeros(M/2,1);
    for iter = 1:M/2
        TempComb1 = CombinationsForNight1(iter,:);
        TempComb2 = CombinationsForNight2(iter,:);
        ScalarProduct(iter) = RatesForNights1(TempComb1)'*RatesForNights2(TempComb2)/norm(RatesForNights1(TempComb1),2)/norm(RatesForNights2(TempComb2));
    end
    SortedScalarProduct = sort(ScalarProduct);
    SortedScalarProduct = SortedScalarProduct(:);
end

function [PermutationBank] = Get_M_Combinations_to_N( M, N )
    PermutationBank = nan(1.5*M, N);
    for iter = 1:1.5*M
        PermutationBank(iter,:) = randperm(N);
    end
    PermutationBank = unique(PermutationBank,'rows');

    if(~(size(PermutationBank,1)<M))
        PermutationBank = PermutationBank(1:M,:);
    end
end

%% Getting the actual distribution:

function pairswiseSP = getPairswiseSP(rateMatrix, nbRecordings)

    AllNightCombinations = nchoosek(1:nbRecordings,2);
    nbNightCombos = size(AllNightCombinations,1);
    pairswiseSP = nan(1,nbNightCombos);
    for iComb = 1:nbNightCombos
        InterVal1 = AllNightCombinations(iComb,1);
        InterVal2 = AllNightCombinations(iComb,2);

        RatePair = rateMatrix(:, [InterVal1, InterVal2]);

         pairswiseSP(iComb) =  RatePair(:,1)*RatePair(:,2)';
    end

end