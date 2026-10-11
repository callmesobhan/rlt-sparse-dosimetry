function T = issueTable(severity, row, message)
%ISSUETABLE Build the standard validation-issue table (Severity, DataRow, Message).
%   Severity is "error" or "warning"; DataRow is the 1-based data row (NaN if
%   the issue is not tied to a row). The column cannot be called Row: MATLAB tables
%   reserve that name for the row dimension. Called with no arguments it returns an empty table.
if nargin == 0
    severity = strings(0, 1);
    row      = zeros(0, 1);
    message  = strings(0, 1);
end
T = table(severity(:), row(:), message(:), ...
    'VariableNames', {'Severity', 'DataRow', 'Message'});
end
