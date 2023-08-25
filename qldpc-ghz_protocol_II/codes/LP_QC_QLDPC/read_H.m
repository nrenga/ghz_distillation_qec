function H = read_H(Hfilename)

f_id = fopen(Hfilename);
line1 = fgetl(f_id);
entries = split(strtrim(line1));
n = str2double(entries{1});
m = str2double(entries{2});

no_of_1s = 0;
for vn = 1:n
    current_line = fgetl(f_id);
    entries = split(strtrim(current_line));
    no_of_1s = no_of_1s + str2double(entries{1});
end

fclose(f_id);

f_id = fopen(Hfilename);
fgetl(f_id);

rows = zeros(no_of_1s,1);
cols = zeros(no_of_1s,1);

ind = 1;
for vn = 1:n
    current_line = fgetl(f_id);
    entries = split(strtrim(current_line));
    vn_deg = str2double(entries{1});
    for i = 1:vn_deg
        cn = str2double(entries{i+1});
        rows(ind) = cn;
        cols(ind) = vn;
        ind = ind + 1;
    end
end

H = sparse(rows,cols,ones(no_of_1s,1));

fclose(f_id);

end