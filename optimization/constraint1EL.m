function loss = constraint1EL(x)
loss = 0;
ALPHA = x(1);
MU = x(2);
KELAS = x(3);
ALPHAVIS = x(4);
MUVIS = x(5);
KVIS = x(6);
RTIME = x(7);
if ALPHA * MU < 0
    loss = 1;
end
if ALPHAVIS * MUVIS < 0
    loss = 1;
end

if KELAS < 0 || KVIS < 0 
    loss = 1;
end
if RTIME < 0 
    loss = 1;
end
end
