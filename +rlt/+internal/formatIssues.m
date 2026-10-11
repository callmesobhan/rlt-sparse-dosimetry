function s = formatIssues(issues, maxLines)
%FORMATISSUES Render validation issues as a readable multi-line string.
%   Shows at most maxLines issues (default 15) and says how many were hidden.
if nargin < 2
    maxLines = 15;
end
n = height(issues);
lines = strings(min(n, maxLines), 1);
for i = 1:numel(lines)
    if isnan(issues.Row(i))
        where = "";
    else
        where = sprintf("row %d: ", issues.Row(i));
    end
    lines(i) = "  - " + where + issues.Message(i);
end
s = strjoin(lines, newline);
if n > maxLines
    s = s + newline + sprintf("  ... and %d more.", n - maxLines);
end
end
