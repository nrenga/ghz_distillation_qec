function print_matrix_with_spaces(M, spaces, long_spaces, hlines)
% Print the rows of matrix M with 'spaces' giving col indices after which a
% space must be printed, 'hlines' giving row indices to be followed with
% a horizontal line

% Author: Narayanan Rengaswamy (July 13, 2021)

if (nargin == 2)
    long_spaces = [];
    hlines = [];
end

if (nargin == 3)
    hlines = [];
end

fprintf('\n');

for i = 1:size(M,1)
    fprintf('\t');
    for j = 1:size(M,2)
        if (M(i,j) >= 0)
            fprintf('%d',M(i,j));
        else
            fprintf('\b%d',M(i,j));
        end
        if (any(spaces==j))
            fprintf(' ');
        end
        if (any(long_spaces==j))
            fprintf('   ');
        end
    end
    fprintf('\n');
    if (any(hlines==i))
        fprintf('\t');
        for j = 1:(size(M,2) + length(spaces) + 3*length(long_spaces))
            fprintf('-');
        end
        fprintf('\n');
    end
end

fprintf('\n');

end