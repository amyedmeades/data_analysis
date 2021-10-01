datalist = dir('test data'); %loads up the list of files in the folder
fields = {'folder' 'date' 'bytes' 'isdir' 'datenum'}; %remove all of the struct fields except the file name
datalist = rmfield(datalist,fields);
datalist = struct2cell(datalist); %convert from structure to cell array
datalist = datalist(3:length(datalist)); %remove the first two blank cells
datalist = datalist'; %transpose
last=length(datalist);

for i= 1:last; %for each data file
    filename = datalist(i); %get the file name from the list
    filename = string(filename); %convert the file name to a string
    filename = char(filename); %convert the string to a character array so we can pick out a section
    data = readmatrix(filename); %read in the data
    w1(i,:)=data(2:257,1); %select w1 data and add to matrix. This could be more automated
    w3(i,:)=data(1,2:365); % select w3 data and add to matrix. This could be more automated
    R(i,:,:)=data(2:257,2:365); %select R data and add to matrix
    time=strfind(filename,'_'); %search filename for underscore to indicate where t2 is
    t2pos=time(1)+1; %calculate the position of t2 in the filename
    t2(i)=str2double(filename(t2pos:t2pos+4)); %extract t2 from the file name and add to matrix
end



