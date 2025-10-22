function [outputArg1,outputArg2] = errorbar_text(X,Y,Xerr,Yerr,textcell,color_order)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if length(Xerr)==1
    Xerr=Xerr*ones(size(X));
end

if length(Yerr)==1
    Yerr=Yerr*ones(size(Y));
end

colors=get(gca,'colororder');
if ~exist('color_order')
    color_order=1:length(textcell);
end

hold on
ii=length(textcell);
errorbar(X(ii+1:end),Y(ii+1:end),Yerr(ii+1:end),Yerr(ii+1:end),...
    Xerr(ii+1:end),Xerr(ii+1:end),'.','MarkerSize',20,'CapSize',0,'color',[.5 .5 .5])
for ii=1:length(textcell)
    errorbar(X(ii),Y(ii),Yerr(ii),Yerr(ii),Xerr(ii),Xerr(ii),'.','MarkerSize',...
        20,'CapSize',0,'color',colors(color_order(ii),:))
end
for ii=1:length(textcell)
    text(X(ii)+.01,Y(ii)+.01,textcell{ii},'color',colors(color_order(ii),:))
end
end

