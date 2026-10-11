function p = provenance()
%PROVENANCE Where and when a result was produced (for reports and audit trails).
p.toolbox   = "rlt-sparse-dosimetry";
p.version   = rlt.version();
p.matlab    = string(version);
p.release   = string(version('-release'));
p.computer  = string(computer);
p.timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end
