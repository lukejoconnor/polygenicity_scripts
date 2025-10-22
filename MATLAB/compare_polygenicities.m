% Compare three polygenicity metrics across a range of possible genetic
% architectures. Reproduces Figure 1 from O'Connor & Sella bioRxiv.

clear;clc
addpath('helpers')

f = @(x)x; finv = f;
g = @log; ginv = @exp;
h = @(x)exp(1/max(x(:)) - 1./x);
hinv = @(y,xmax)-1./(-1/xmax + log(y));

Pi = {@(x)compute_true_polygenicity(x,x,g,ginv);
    @(x)compute_true_polygenicity(x,x,f,finv);
    @(x)compute_true_polygenicity(x,x,h,hinv);
    @(x)sum(x>0)};

% small effects to zero with fixed h2
large_effect = 0.01;
n_large = 10;
n_small_array = 10*2.^(0:10);
for ii=1:length(n_small_array)
    n_small = n_small_array(ii);
    small_effect = large_effect * n_large / n_small;
    x = [large_effect * ones(n_large,1); small_effect * ones(n_small,1)];
    for p = 1:length(Pi)
        poly(ii,p) = Pi{p}(x);
    end
    
end

figure;subplot(1,4,2)
h = plot(n_small_array, poly);
labels = {'{\it P_G} (entropy)',...
    '{\it P_H} (effective)',...
    '{\it P_M} (softmax)',...
    '{\it P_A} (nonzero effects)'};
order = [4,1,2,3];
legend(h(order),labels(order))
legend boxoff
set(gca,'yscale','log','xscale','log')
set(gca,'ytick',[10, 100,1000])
xlabel('Number of small effects (n_s)')
ylabel('Polygenicity')
ylim([10 5000])
title('A')
box off

% fraction of small effect h2 to 1
clear poly
large_effect = 0.01;
small_effect = 0.0001;
n_large_array = 0:20;
n_small_array = round((20 - n_large_array) * large_effect / small_effect);
for ii=1:length(n_small_array)
    n_small = n_small_array(ii);
    n_large = n_large_array(ii);
    x = [large_effect * ones(n_large,1); small_effect * ones(n_small,1)];
    for p = 1:length(Pi)
        poly(ii,p) = Pi{p}(x);
    end
    
end
fraction_small = n_small_array * small_effect;
fraction_small = fraction_small ./ ...
    (fraction_small + n_large_array * large_effect);

subplot(1,4,3)
plot(fraction_small, poly)
set(gca,'yscale','log')
set(gca,'ytick',[10, 100,1000])
xlabel('Contribution of small effects (n_sp_s / (n_sp_s + n_lp_l))')
ylabel('Polygenicity')
ylim([10 5000])
title('B')
box off

clear poly
n_small = 1e3;
n_large = 10;
large_effect = .01;
small_effect_array = large_effect * (0:10:1000)/1000;
for ii=1:length(small_effect_array)
    small_effect = small_effect_array(ii);
    x = [large_effect * ones(n_large,1); small_effect * ones(n_small,1)];
    for p = 1:length(Pi)
        poly(ii,p) = Pi{p}(x);
    end
end
ratio_small = small_effect_array / large_effect;

subplot(1,4,4)
plot(ratio_small', poly)
set(gca,'yscale','log','xscale','linear')
set(gca,'ytick',[10, 100, 1000])
xlabel('Contribution per small effect (p_s/p_l)')
ylabel('Polygenicity')
ylim([10 5000])
title('C')
box off
