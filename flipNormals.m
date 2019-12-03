    function flipNormals(filePattern,outputSuffix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   By Daniel Tward, Aug. 26, 2013
%
%   inputs: 
%       filePattern     -   a pattern for the dir command.  Can be a
%                           filename, or something like "*.byu".
%       outputSuffix    -   OPTIONAL.  A string to append to the filename.
%                           The default is "flipped".
%   requirements:
%       loadbyu_local and writebyu must be in your current directory or
%       matlab path.
%
%

if nargin == 1
    outputSuffix = 'flipped';
end

if nargin == 0
    error('Must specify at least one argument, file pattern string');
end


% get the files to flip
files = dir(filePattern);

% loop through the files
for i = 1 : length(files)
    % get the current file
    byufile = files(i).name;
    
    % load the surface.  Make sure this is in 
    [nvertex, ntris, nconns, triloc, tris] = loadbyu_local(byufile);
    
    % flip the normals
    tris = [tris(:,2) tris(:,1) tris(:,3)];
    
    % get the new filename
    if outputSuffix(1) == '_';
        separator = '';
    else 
        separator = '_';
    end
    outfile = [byufile(1:end-4) separator outputSuffix '.byu'];
    
    
    % write it
    writebyu(nvertex, ntris, nconns, triloc, tris, outfile);
end