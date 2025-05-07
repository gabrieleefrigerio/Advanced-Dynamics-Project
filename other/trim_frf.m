clear
clc

load("data\Data.mat");

fcut = 1800;
index = find(freq<fcut);
frf = frf(index, :);
cohe = cohe(index, :);
freq = freq(index, :);

save("data\Data_trim.mat", "frf", "cohe", "freq", "fcut");
