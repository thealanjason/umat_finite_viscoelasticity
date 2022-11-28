function loss = constraint3EL(x)
loss = 0;
ALPHA = x(1);
MU = x(2);
KELAS = x(3);
ALPHAVIS = x(4);
MUVIS = x(5);
KVIS = x(6);
RTIME = x(7);
ALPHAVIS_2 = ALPHAVIS;
MUVIS_2 = MUVIS;
KVIS_2 = KVIS;
RTIME_2 = x(8);
ALPHAVIS_3 = ALPHAVIS;
MUVIS_3 = MUVIS;
KVIS_3 = KVIS;
RTIME_3 = x(9);
if ALPHA * MU < 0
    loss = 1;
end
if ALPHAVIS * MUVIS < 0
    loss = 1;
end
if ALPHAVIS_2 * MUVIS_2 < 0
    loss = 1;
end
if ALPHAVIS_3 * MUVIS_3 < 0
    loss = 1;
end

if KELAS < 0 || KVIS < 0 || KVIS_2 < 0 || KVIS_3 < 0
    loss = 1;
end
if RTIME < 0 || RTIME_2 < 0 || RTIME_3 < 0
    loss = 1;
end
end
