function data = ReadFile( file, num_elem )
fileID = fopen(file,'r');
if fileID == -1
    error('Could not open file')
end
tline = fgetl(fileID); %Not using descriptions

line_format = repmat('%f,',1,num_elem);
line_format = line_format(1:end-1);

j = 1;
%n = i_end - i_start;
f = [0,0];
while 1
    tline = fgetl(fileID);
    if ~ischar(tline), break, end
    %disp(tline)
    
    f = sscanf(tline, line_format);
    
    data(j,1:num_elem) = f;
    j = j + 1;
    
end

fclose(fileID);

end