function databrushlabels =  regionNamestoNo(databrushreg)
if ischar(databrushreg)
    databrushreg = {databrushreg};
%     error('Should not input char array to databrushreg')
end

databrushlabels = zeros(length(databrushreg),1);

indices = find(strcmp(databrushreg, 'ManRO'));
databrushlabels(indices) = 1*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'ManRB'));
databrushlabels(indices) = 2*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'ManRL'));
databrushlabels(indices) = 3*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'ManAB'));
databrushlabels(indices) = 4*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'ManAL'));
databrushlabels(indices) = 5*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'ManLO'));
databrushlabels(indices) = 6*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'ManLB'));
databrushlabels(indices) = 7*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'ManLL'));
databrushlabels(indices) = 8*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxRO'));
databrushlabels(indices) = 9*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxRB'));
databrushlabels(indices) = 10*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxRL'));
databrushlabels(indices) = 11*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxAB'));
databrushlabels(indices) = 12*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxAL'));
databrushlabels(indices) = 13*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxLO'));
databrushlabels(indices) = 14*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxLB'));
databrushlabels(indices) = 15*ones(length(indices),1);
indices = find(strcmp(databrushreg, 'MaxLL'));
databrushlabels(indices) = 16*ones(length(indices),1);


end