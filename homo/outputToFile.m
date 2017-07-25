function outputToFile(coeffest, outfname)

fout = fopen(outfname, 'w');
for i=1:size(coeffest,1)
    arr = coeffest(i,:);
    fprintf(fout, '%s %s %s\n', arr{:});
end;
fclose(fout);