function [] = matrix2csv(matrix, filename, colnames, rownames, delimiter)
%MATRIX2CSV Write matrix to CSV file.
%   MATRIX2CSV(MATRIX, FILENAME) writes MATRIX to a CSV file.
%
%   MATRIX2CSV(MATRIX, FILENAME, COLNAMES) includes column headers.
%
%   MATRIX2CSV(MATRIX, FILENAME, COLNAMES, ROWNAMES) includes row names.
%
%   MATRIX2CSV(..., DELIMITER) uses the specified delimiter (default: ',').

% parse inputs
if nargin < 3
    colnames = {}; rownames = {}; delimiter = ',';
elseif nargin < 4
    rownames = {}; delimiter = ',';
elseif nargin < 5
    delimiter = ',';
end

colnames_flag = ~isempty(colnames);
rownames_flag = ~isempty(rownames);
if colnames_flag && rownames_flag
    table = array2table(matrix, "RowNames", rownames, "VariableNames", colnames);
elseif colnames_flag && ~rownames_flag
    table = array2table(matrix, "VariableNames", colnames);
elseif ~colnames_flag && rownames_flag
    table = array2table(matrix, "RowNames", rownames);
else
    table = array2table(matrix);
end

writetable(table, filename, "Delimiter", delimiter, "WriteVariableNames", colnames_flag);
end
