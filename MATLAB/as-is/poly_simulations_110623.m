% This simulation script is provided as-is for completeness.
% Running it will require manual path manipulation as well as
% installing SuiteSparse, cloning the LDGMs repo, and 
% downloading UK Biobank LDGM precision matrices.

clear;clc
if ~exist('P')
    indexFn = @(x,j)x(j);
    id = @(x)x;

    addpath(genpath('/Users/loconnor/Library/CloudStorage/Dropbox/GitHub/ldgm/MATLAB/'))
    addpath(genpath('/Users/loconnor/SuiteSparse/MATLAB_Tools'))

    chr = 1;
    populations = 'EUR';
    ldgms_dir = '~/Documents/data/precision_ukbb/';
    % note that SNP lists are 1kg, but precision matrices are UKBB
    snplist_dir = ['~/Documents/data/ldgm/1kg_chr', num2str(chr), '_'];
    [P, snplists, AF] = loadLDGMs(ldgms_dir,'popnNames','EUR','snplistpath',snplist_dir);
    noBlocks=length(P);

    % set AF to zero for SNPs not in LDGM
    for block = 1:noBlocks
        idx = any(P{block});
        AF{block}(~idx) = 0;
    end

    % Correlation matrices
    R = cell(size(P));
    for block = 1:noBlocks
        idx = any(P{block});
        R{block} = full(inv(P{block}(idx,idx)));
    end
    
    % Annotation matrices to allow multiple SNPs per index
    whichIndicesAnnot = cell(noBlocks,1);
    annot = whichIndicesAnnot;
    for block = 1:noBlocks
        idx = find(any(P{block}));
        x = snplists{block}.index+1;
        whichIndicesAnnot{block} = x(ismember(x,idx));
        af = AF{block}(whichIndicesAnnot{block});
        annot{block} = min(af,1-af) > 0.05;
    end

end

T = table();
T.description{1} = 'Default';

% Simulation parameters
T.sampleSize = 1e5;
T.heritability = 0.1;
T.alpha = -1;

sparsity = [200 500 1000 2000 5000 10000];
noRows = length(sparsity);
T = repmat(T,noRows,1);
for row = 1:noRows
    T.componentWeight(row,:) = {[1 0.2] / sparsity(row) };
    T.componentVariance(row,:) = {1./T.componentWeight{row,:}};
end

noReplicates = 100;
T = repmat(T,noReplicates,1);

poly_functions = {@log, id, @(x)exp(1/max(x(:)) - 1./x)};
inverse_functions = {@exp, id, @(y,xmax)-1./(-1./xmax + log(y))};

for counter = 1:noRows*noReplicates
    disp(counter)
    % simulate summary statistics
    tic;
    rng(ceil(counter/noRows),'twister')
    [sumstats, whichIndicesSumstats, ~, true_beta_perSD, ~, perSNPh2] = ...
        simulateSumstats(T.sampleSize(counter),...
        'precisionMatrices',P,...
        'whichIndicesAnnot',whichIndicesAnnot,...
        'annotations',annot,...
        'linkFn',@(x)1,...
        'alleleFrequency', AF,...
        'heritability', T.heritability(counter),...
        'componentWeight', T.componentWeight{counter},...
        'componentVariance', T.componentVariance{counter},...
        'alphaParam', T.alpha(counter),...
        'componentRandomSeed',1,...
        'returnSigmasq',true);
    timeSimulate = toc;

    if counter <= noRows
        marginal_variance = cell(noBlocks,1);
        for block = 1:noBlocks
            marginal_variance{block} = sum(R{block}.^2 .* perSNPh2{block})';
        end

        perSNPh2Cat = vertcat(perSNPh2{:});
        marginal_variance = vertcat(marginal_variance{:});
        for ii = 1:length(poly_functions)
            T.true_polygenicity(counter,ii) = ...
                compute_true_polygenicity(perSNPh2Cat,marginal_variance,...
                poly_functions{ii},inverse_functions{ii});
        end
    else
        T.true_polygenicity(counter,:) = T.true_polygenicity(mod(counter-1,noRows)+1,:);
    end

    % Z scores
    T.Z{counter} = cellfun(@(T)T.Z_deriv_allele,sumstats,'UniformOutput',false);
    T.whichIndicesSumstats{counter} = whichIndicesSumstats;
    disp(T.true_polygenicity(counter,:))
end


