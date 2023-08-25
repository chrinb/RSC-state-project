function plot_swr(sData)

time = (1:length(sData.ephysdata.rippleSnips(149).lfp))/2500;

for i = 1:length(sData.ephysdata.absRipIdx)

    figure,
    plot(time, sData.ephysdata.rippleSnips(i).lfp)
    set(gca,'xlim',[time(1) time(end)])
    prompt = sprintf('Next SWR? ');
    x       = input(prompt,'s');
    if isempty(x)
        close 
        clc
    end
end
