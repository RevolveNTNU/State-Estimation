function ch = import_channel_analyze(filename)

%filename = 'C:\Users\tonja\OneDrive\Skrivebord\CSV\FSG - AutoX run 1 and 2\vcu.INS.ax.csv';
[fid, msg] = fopen(filename, 'rt');
if fid < 0
    error('Could not open file "filename" because: "%s"', filename, msg);
end
ch = {};
while true
    thisline = fgetl(fid);
    if ~ischar(thisline); break; end   %normal: indicates end of file
    if isempty(thisline); continue; end  %empty line, including possibly just before end of file
    fields = regexp(thisline, ',', 'split');
    if length(fields) >= 2
       ch{end+1} = fields{2}; 
    end
end
fclose(fid);
temp = str2double(ch);
if ~any(isnan(temp))
    ch = temp;    %convert to double if ALL entries can be converted
end

end

