function L = read_L(Lfilename)

f_id = fopen(Lfilename);
line1 = fgetl(f_id);
entries = split(strtrim(line1));
k = str2double(entries{1});
n = str2double(entries{2});

line2 = fgetl(f_id);
entries = split(strtrim(line2));
row_wt_max = str2double(entries{1});
col_wt_max = str2double(entries{2});

line3 = fgetl(f_id);
entries = split(strtrim(line3));
row_wts = cellfun(@str2double,entries);

line4 = fgetl(f_id);
entries = split(strtrim(line4));
col_wts = cellfun(@str2double,entries);

if (sum(row_wts) ~= sum(col_wts))
    fprintf('\n Sum of row weights NOT EQUAL TO sum of column weights! \n');
    L = [];
    return;
end

no_of_1s = sum(row_wts);
rows = zeros(no_of_1s,1);
cols = zeros(no_of_1s,1);

ind = 0;
for q = 1:k
    current_line = fgetl(f_id);
    entries = split(strtrim(current_line));
    row_supp = cellfun(@str2double,entries);
    rows((ind+1):(ind+row_wts(q))) = q * ones(row_wts(q),1);
    cols((ind+1):(ind+row_wts(q))) = row_supp;
    ind = ind + row_wts(q);
end

L = sparse(rows,cols,ones(no_of_1s,1));

fclose(f_id);

end