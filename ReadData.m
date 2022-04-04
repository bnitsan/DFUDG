
function [tab] = ReadData()

tab1=readSingleFile('../gal1_cat_2reff.dat');
tab2=readSingleFile('../gal2_cat_2reff.dat');

tab = [tab1' tab2']';
%save('dat.mat','tab')
end

function tab=readSingleFile(name)

tab0=readtable(name);

xid  = tab0.Var1;
id   = tab0.id;
RA   = tab0.ra;
Dec  = tab0.dec;  
f606tot = tab0.f606w_total_mag;  
f475tot = tab0.f475w_total_mag;  
vel     = tab0.id*0;
vel_err = tab0.id*0;    
specconf= tab0.spec_conf;  
f606err = tab0.id*0;
f475err = tab0.id*0;

tab = [xid id RA Dec f606tot f475tot vel vel_err specconf f606err f475err];
end
