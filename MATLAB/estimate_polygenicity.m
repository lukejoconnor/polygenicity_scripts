% Aug 29 2024
% Estimate polygenicity of 36 traits using multiple related definitions

clear;clc
matfiles_path = '../../FMR/matfiles/';
fmr_path = '../../FMR/';
addpath(genpath([fmr_path, 'MATLAB']))
addpath('helpers')
table_save_path = '../polygenicity_estimates.txt';

% Load data
FMR1 = load([matfiles_path, 'FMRestimates_32traits.mat'],...
    'ss_est','ww_jk','LD4Mout','mm','traits',...
    'h2gs_obs','warningflag','numgs_obs');

FMR2 = load([matfiles_path, 'FMR_lipids.mat'],...
    'ss_est','ww_jk','LD4Mout','mm','traits',...
    'h2gs_obs','warningflag','numgs_obs');

load([matfiles_path, 'name_abbreviations.mat'], 'abbrevs', 'traits')
load([matfiles_path, '32wellpowered_sumstats.mat'], 'sumstats')

% Process input data
for tt = 1:length(sumstats)
    maxchisq(tt) = max(sumstats(tt).chisq_N145k);
    h2(tt) = mean(FMR1.LD4Mout(tt).cov) * FMR1.mm;
    h2_over_max(tt) = h2(tt) / maxchisq(tt);
end
mm = FMR1.mm;
incl = [1:3 5]; % new traits: exclude skin color, which is very noisy
ss_est = [FMR1.ss_est; FMR2.ss_est(incl,:)];
ww_jk = [FMR1.ww_jk; FMR2.ww_jk(incl,:,:)];
LD4Mout = [FMR1.LD4Mout'; FMR2.LD4Mout(incl)'];
h2gs_obs = [FMR1.h2gs_obs; FMR2.h2gs_obs(incl,:)];
numgs_obs = [FMR1.numgs_obs; FMR2.numgs_obs(incl,:)];
nt = length(traits);
clear FMR1 FMR2 sumstats

%% Estimate polygenicity

% Polygenicity functions
logsumexp = @(a,b)max(a,b) + log(1+exp(-abs(b-a)));
logdifexp= @(a,b)a + log(1-exp(b-a));
poly_functions = {@log, @(x)x, @(x)max(1e-256,exp(-1./x))};
inverse_functions = {@exp, @(x)x, @(x)-1./log(x)};

% Estimates
poly_jk = zeros(100,length(poly_functions));
for tt=1:nt
    disp(traits{tt})
    ww = reshape(ww_jk(tt,:,:),13,100)';
    Nh2 = mm * mean(LD4Mout(tt).cov);
    
    for ii = 1:length(poly_functions)
        poly_jk(:,ii) = compute_polygenicity(ss_est(tt,:)/Nh2,ww,poly_functions{ii},inverse_functions{ii});
    end

    poly_est(tt,:) = mean(log10(poly_jk));
    poly_err(tt,:) = std(log10(poly_jk)) * sqrt(102);

end

%% Plotting

special_plot_pheno=[1 6 14 28 34 36];
special_pheno_names={'BMI','Height','Schizophrenia','IBD','LDL','Hair color'};

% Figure 2B
figure;subplot(3,2,2)
plot([0 5],[0 5], 'color', [.8 .8 .8]); hold on
order = [special_plot_pheno, setdiff(1:nt,special_plot_pheno)];
errorbar_text(poly_est(order,1),poly_est(order,2),poly_err(order,1),poly_err(order,2),special_pheno_names)
xlabel('Entropy polygenicity (SE)')
ylabel('Effective polygenicity (SE)')
set(gca,'YTick',1:5,'YTickLabel',{'10','100','1000','10^4','10^5'})
set(gca,'XTick',1:5,'XTickLabel',{'10','100','1000','10^4','10^5'})
ylim([0.5 5])
xlim([2 5.5])
title('B')
box off

% Figure 2C
subplot(3,2,4)
plot([0 5],[0 5], 'color', [.8 .8 .8]); hold on
errorbar_text(poly_est(order,1),poly_est(order,3),poly_err(order,1),poly_err(order,3),special_pheno_names)
xlabel('Entropy polygenicity (SE)')
ylabel('Softmax polygenicity (SE)')
set(gca,'YTick',1:5,'YTickLabel',{'10','100','1000','10^4','10^5'})
set(gca,'XTick',1:5,'XTickLabel',{'10','100','1000','10^4','10^5'})
ylim([0.5 5])
xlim([2 5.5])
title('C')
box off

% Figure 2E
subplot(3,2,6)
errorbar_text(poly_est(order,1), numgs_obs(order,1), poly_err(order,1), zeros(nt,1), special_pheno_names)
set(gca,'XTick',1:5,'XTickLabel',{'10','100','1000','10^4','10^5'})
xlabel('Entropy polygenicity (SE)')
xlim([2 5.5])
ylabel('No. sig loci')
title('E')
box off

% Figure 2F
subplot(3,2,5)
match_targetsize_poly
box off
title('D')

% Figure 2A
colors = colororder;
[~,order] = sort(poly_est(:,1));
subplot(3,2,1); hold on; clear h
for ii=1:3
    h(ii) = errorbar(poly_est(order,ii),poly_err(order,ii),'o','markerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:));
end
set(gca,'XTick',1:nt,'XTickLabel',abbrevs(order))
set(gca,'YTick',0:5,'YTickLabel',{'1','10','100','1000','10^4','10^5'})
ylim([0 5])
view(90,-90)
labels = {'\Pi_G (entropy)', '\Pi_H (effective)', '\Pi_G (softmax)'};
legend(h, labels)
legend boxoff
ylabel('Polygenicity (SE)')
xlim([0.5 nt + .5])
title('A')
box off

% Supp figure with maximum contribution per trait
figure; hold on
plot([1,3],[1,3], 'black')
scatter(poly_est(1:length(h2_over_max),3), log10(h2_over_max), ...
    'o','markerFaceColor',colors(3,:),'markerEdgeColor',colors(3,:));
set(gca,'XTick',1:3,'XTickLabel',{'10','100','1000'})
set(gca,'YTick',1:3,'YTickLabel',{'10','100','1000'})
xlabel('Entropy polygenicity (SE)')
xlim([1,3]); ylim([1,3])
xlabel('Softmax polygenicity')
ylabel('h^2 / max(p)')

% Results table
entropy_polygenicity = poly_est(:,1);
effective_polygenicity = poly_est(:,2);
softmax_polygenicity = poly_est(:,3);
entropy_polygenicity_se = poly_err(:,1);
effective_polygenicity_se = poly_err(:,2);
softmax_polygenicity_se = poly_err(:,3);
num_gws_loci = numgs_obs(:,1);
target_size_simons = zeros(nt,1);
target_size_simons_lb = zeros(nt,1);
target_size_simons_ub = zeros(nt,1);
target_size_simons(matches(:,1)) = target_size(matches(:,2));
target_size_simons_lb(matches(:,1)) = target_size_bounds(matches(:,2),1);
target_size_simons_ub(matches(:,1)) = target_size_bounds(matches(:,2),2);
h2_over_max(length(h2_over_max):length(entropy_polygenicity)) = 0;
h2_over_max = reshape(h2_over_max, length(h2_over_max), 1);

T = table(traits, entropy_polygenicity, entropy_polygenicity_se,...
    effective_polygenicity, effective_polygenicity_se,...
    softmax_polygenicity, softmax_polygenicity_se, ...
    num_gws_loci,...
    h2_over_max,...
    target_size_simons, target_size_simons_lb, target_size_simons_ub);

writetable(T, table_save_path)

