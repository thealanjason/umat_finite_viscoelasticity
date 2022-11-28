function loss = constraint2EL(x)
loss = 0;
ALPHA = x(1);
MU = x(2);
KELAS = x(3);
ALPHAVIS = x(4);
MUVIS = x(5);
KVIS = x(6);
RTIME = x(7);
ALPHAVIS_2 = x(8);
MUVIS_2 = x(9);
KVIS_2 = x(10);
RTIME_2 = x(11);
if ALPHA * MU < 0
    loss = 1;
end
if ALPHAVIS * MUVIS < 0
    loss = 1;
end
if ALPHAVIS_2 * MUVIS_2 < 0
    loss = 1;
end

if KELAS < 0 || KVIS < 0 || KVIS_2 < 0
    loss = 1;
end
if RTIME < 0 || RTIME_2 < 0
    loss = 1;
end
end
