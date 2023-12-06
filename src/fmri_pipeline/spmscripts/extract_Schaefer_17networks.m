function extract_Schaefer_17networks(FCfile)

% Date created: 2022-01-20, Author: Hweeling Lee
% this script creates the mean FC for each of the 17 networks. Negative
% correlations are set to 0, before calculating mean FC.

%FCfile='../FC_Schaefer_17Networks_200p_fullcorr.csv';

fs = filesep;

[path, fn, ext] = fileparts(FCfile);

RSNetworkNames = {'ContA', 'ContB', 'ContC', 'DefaultA', 'DefaultB', 'DefaultC', 'DAttnA', 'DAttnB', 'LimbicA', 'LimbicB', 'SVAttnA', 'SVAttnB', 'SomMotA', 'SomMotB', 'TempPar', 'VisCent', 'VisPeri'};
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


end

