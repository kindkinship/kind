function [map,ped,frq,samplePairsKnown,coeffKnown] = readInputFiles(mapfile, pedfile, frqfile, knownrel)

map = {};
ped = {};
frq = [];
samplePairsKnown = {};
coeffKnown = [];

% read the map file
fid = fopen(mapfile,'r');
C = textscan(fid, '%s %s %s %s');
fclose(fid);
map = [C{1} C{2} C{3} C{4}];

% check whether there are non-autosomal chromosomes in the data
uniqueChrnos = unique(map{1});
if ~isempty(find(strcmpi('x',uniqueChrnos))) || ~isempty(find(strcmpi('y',uniqueChrnos)))
    error('Only autosomal chromosomes are supported.');
end;

% read the ped file
formatStr = ['%s %s %s %s %s %s' repmat(' %s',1,2*size(map,1))];
fid = fopen(pedfile,'r');
C = textscan(fid, formatStr);
fclose(fid);

ped = {};
for i=1:size(C,2)
    ped = [ped C{i}];
end;    

% read the allele frequency file
fid = fopen(frqfile,'r');
C = textscan(fid, '%s %s %s %s %s %s');
fclose(fid);
frq = str2double(C{5});

% read the known relationship file
fid = fopen(knownrel,'r');
C = textscan(fid, '%s %s %f');
fclose(fid);
samplePairsKnown = [C{1} C{2}];
coeffKnown = C{3};
