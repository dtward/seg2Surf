function [f,v] = readbyu_local(filename)

% daniel wrote this

% open file
fid = fopen(filename);

% get the first line
line1 = fgetl(fid);
[T,R] = strtok(line1); % ignore
[T,R] = strtok(R);
nv = str2num(T);
[T,R] = strtok(R);
nf = str2num(T);

% ignore the second line
line2 = fgetl(fid);

% read vertices
v = fscanf(fid,'%f %f %f\n',[3 nv])';

% read faces
f = abs(fscanf(fid,'%d %d %d\n',[3 nf])');

% close file
fclose(fid);

