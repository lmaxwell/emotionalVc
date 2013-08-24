function frequencyAnchors = setBlankFrequencyAnchors(label)

upperLimitOfFrequencyAnchors = 6;
frequencyAnchors.counts = zeros(size(label.segment,1),1);
frequencyAnchors.frequency = ...
    zeros(size(label.segment,1),upperLimitOfFrequencyAnchors);
for ii = 1:size(label.segment,1)
    for jj = 1:upperLimitOfFrequencyAnchors
        markA = -500+jj*1000-50;
        frequencyAnchors.frequency(ii,jj) = markA;
    end;
end;

