function T = readCsvTable(file, textVars, numericVars)
%READCSVTABLE Read a CSV with declared column types instead of guessed ones.
%   T = rlt.internal.readCsvTable(file, textVars, numericVars)
%
%   readtable infers each column's type from its contents, which goes wrong on
%   exactly the data this toolbox imports: a column mixing "yes" and 0 becomes
%   numeric with NaN for "yes"; an ID column of digits becomes a number; a header
%   row can be mistaken for data, and the delimiter itself can be guessed wrongly
%   when a quoted text field is long. Here the delimiter is always a comma, line 1 is
%   always the header, columns named
%   in textVars are read as string, and columns named in numericVars as double.
%   Names not present in the file are ignored (validators report missing columns).
arguments
    file (1,1) string
    textVars (1,:) string
    numericVars (1,:) string
end
opts = detectImportOptions(file, 'FileType', 'delimitedtext', 'Delimiter', ',');
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
present = string(opts.VariableNames);
for v = textVars(ismember(textVars, present))
    opts = setvartype(opts, char(v), 'string');
end
for v = numericVars(ismember(numericVars, present))
    opts = setvartype(opts, char(v), 'double');
end
T = readtable(file, opts);
end
