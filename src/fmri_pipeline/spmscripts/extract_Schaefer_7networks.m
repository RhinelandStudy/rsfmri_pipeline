function extract_Schaefer_7networks(FCfile)

%this script creates the mean FC for each of the 7 networks
%Negative correlations are set to 0, before calculating mean FC.
%FCfile='../FC_Schaefer_7Networks_200p_fullcorr.csv';

fs = filesep;

[path, fn, ext] = fileparts(FCfile);

RSNetworkNames = {'Cont', 'Default', 'DAttn', 'Limbic', 'SVAttn', 'SomMot', 'Vis'};
opts = detectImportOptions(FCfile);
clear tmp1 tmp2 tmp3 RS_networks TID;
tmp1 = regexp(opts.VariableNames', {'_'}, 'split');
for i = 1:size(tmp1,1)
    tmp2(i,:) = [i tmp1{i}(:,2:3)];
end

[tmp4 TID] = findgroups(tmp2(:,2));

for i = 1:size(TID,1)
    clear tmp5;
    tmp5 = find(tmp4 == i);
    RS_networks(i,:) = [tmp5(1,1) tmp5(end,1)];
end

FCMatrix = dlmread(FCfile, ',', 1,0);
FCMatrix = atanh(FCMatrix);
FCMatrix(FCMatrix < 0) = 0;

count = 1;
for i = 1:size(RS_networks,1)
    clear tmp5;
    tmp5 = FCMatrix(RS_networks(i,1):RS_networks(i,2),RS_networks(i,1):RS_networks(i,2));
    FCM(count,1) = mean(tmp5(find(triu(tmp5,1))));
    FCnames{count,1} = [RSNetworkNames{i}, '_', RSNetworkNames{i}];
    count = count + 1;
end

for i = 1:size(RS_networks,1)
    for ii = 1:size(RS_networks,1)
        if ii < i
            clear tmp5;
            tmp5 = FCMatrix(RS_networks(ii,1):RS_networks(ii,2),RS_networks(i,1):RS_networks(i,2));
            FCM(count,1) = mean(tmp5, 'all');
            FCnames{count,1} = [RSNetworkNames{ii}, '_', RSNetworkNames{i}];
            count = count + 1;
        end
    end
end

t2 = [cell2table(FCnames, 'VariableNames', {'FCName'}), table(FCM, 'VariableNames', {'meanFC'})];
disp('writing output csv...');
writetable(t2, strcat(pwd, fs, [fn '_means.csv']));

t2js = jsonencode(t2);
fid=fopen([pwd  fs fn '_means.json'],'w');
disp('writing output json...');
fprintf(fid,t2js);
fclose(fid);

% system segregation
for i = 1:size(RS_networks,1)
    clear within_network_idx between_networks_idx;
    within_network_idx = [RS_networks(i,1):RS_networks(i,2)];
    between_networks_idx = [RS_networks(1,1):RS_networks(7,2)];
    between_networks_idx(within_network_idx) = [];
    
    clear tmp6 tmp7;
    tmp6 = FCMatrix(within_network_idx, within_network_idx);
    tmp7 = FCMatrix(within_network_idx, between_networks_idx);
    
    clear within_network between_network;
    within_network = mean(tmp6(find(triu(tmp6,1))));
    between_network = mean(tmp7, 'all');
    
    system_seg(i,1) = (within_network - between_network) / within_network;
    segnames{i,1} = ['sysseg_', RSNetworkNames{i}];
    
end

t3 = [cell2table(segnames), array2table(system_seg, 'VariableNames', {'system_seg'})];
writetable(t3, strcat(pwd, fs, 'systemseg_7networks.csv'));


end




