function T = issueTable(severity, row, message)
%ISSUETABLE Build the standard validation-issue table (Severity, Row, Message).
%   Severity is "error" or "warning"; Row is the 1-based data row (NaN if the
%   issue is not tied to a row). Called with no arguments it returns an empty table.
if nargin == 0
    severity = strings(0, 1);
    row      = zeros(0, 1);
    message  = strings(0, 1);
end
T = table(severity(:), row(:), message(:), ...
    'VariableNames', {'Severity', 'Row', 'Message'});
end
