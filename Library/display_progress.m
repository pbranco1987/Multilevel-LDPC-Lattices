function display_progress(j,n,step)
if ~exist('step','var')
    step = 10;
end
if n < round(100/step)
    step = 100/n;
end
percent = find(j == ceil(n*[1:round(100/step)]/round(100/step)),1,'last');
if percent
    fprintf([repmat('\b',1,4*(percent>1)) '%2d0%%' repmat('\b',1,4*(percent==round(100/step)))],percent)
end
