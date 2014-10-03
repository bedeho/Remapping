v=[29	7 31 6 34	52 36	25 39	23 41	54 44	65 46	34 48	35 51	64 53	26 56	71 58	38 61	34 63	50 66	81]

for i=v,
    
    x = i-46;
    plot([x x], [0 1]);
    hold on;
end
ylim([0 1]);
xlim([-45 45]);
xlabel('Receptive Field Location');