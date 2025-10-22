function show_eeg(values,fs,labels)


figure
set(gcf,'position',[10 10 1200 800])
dur = size(values,1)/fs;
nchs = length(labels);

offset = 0;
ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for ich = 1:nchs
    
    if all(isnan(values(:,ich)))
        continue
    end
    
    plot(linspace(0,dur,size(values,1)),values(:,ich) - offset,'k');
    hold on
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(values(:,ich));
    
    text(dur+0.05,ch_bl(ich),sprintf('%s',labels{ich}),'fontsize',20)
    
    %{
    if ich < nchs
        if isnan((min(values(:,ich)) - max(values(:,ich+1))))
        else
            offset = offset - (min(values(:,ich)) - max(values(:,ich+1)));
        end    
        
    end
    %}

    if ich < nchs
        jch = ich+1;
        while jch < nchs
            gap = min(values(:,ich)) - max(values(:,jch));
            if isnan(gap)
                jch = jch + 1;
            else
                offset = offset - gap;
                break
            end
        end
    end
    
end

xlabel('Time (seconds)')
set(gca,'fontsize',20)

end