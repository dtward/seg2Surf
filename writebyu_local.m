%function writebyu(nvertex, ntris, nconns, triloc, tris, byufile);
%
% input:
%          nvertex   number of vertices
%          ntris     number of triangles
%          nconns    number of connectivity entries
%          triloc    matrix (#vertices,3) coordinates
%          tris      matric (#triangles,3) verticex in each triangle
%          byufile   byufile name
%          
% this version from daniel, november 5, 2015
function writebyu_local(f,v, byufile);
triloc = v;
tris = f;
nvertex = size(v,1);
ntris = size(f,1);
nconns = ntris*3;

fid = fopen(byufile,'w');
fprintf(fid, '%d ', 1);
fprintf(fid, '%d ', nvertex);
fprintf(fid, '%d ', ntris);
fprintf(fid, '%d\n', nconns);
fprintf(fid, '%d ', 1);
fprintf(fid, '%d\n', ntris);

for i=1:nvertex,
  fprintf(fid, '%f ', triloc(i,1:2));
  fprintf(fid, '%f\n', triloc(i,3));
end

for i=1:ntris,
  fprintf(fid, '%d ', tris(i,1:2));
  fprintf(fid, '%d\n',-tris(i,3));
end


% fclose all;
fclose(fid);% daniel changed to fclose(fid) January 2012



