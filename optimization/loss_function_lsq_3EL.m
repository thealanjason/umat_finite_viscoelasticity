function loss = loss_function_lsq_3EL(x)
    global xsol
    xsol = x;
    MU  = x(1);
    ALPHA = x(2);
    KELAS = x(3);
    MUVIS = x(4);
    ALPHAVIS = x(5);
    KVIS = x(6);
    RTIME = x(7);
    ETADEV = RTIME * MUVIS * ALPHAVIS;
    ETAVOL = RTIME * KVIS;
    MUVIS_2 = x(8);
    ALPHAVIS_2 = x(9);
    KVIS_2 = x(10);
    RTIME_2 = x(11);
    ETADEV_2 = RTIME_2 * MUVIS_2 * ALPHAVIS_2;
    ETAVOL_2 = RTIME_2 * KVIS_2;
    MUVIS_3 = x(12);
    ALPHAVIS_3 = x(13);
    KVIS_3 = x(14);
    RTIME_3 = x(15);
    ETADEV_3 = RTIME_3 * MUVIS_3 * ALPHAVIS_3;
    ETAVOL_3 = RTIME_3 * KVIS_3;
    

    global pt1 pt06 pt03 pt01 pt003
    lenpt1 = size(pt1);
    lenpt1 = lenpt1(1);
    lenpt06 = size(pt06);
    lenpt06 = lenpt06(1);
    lenpt03 = size(pt03);
    lenpt03 = lenpt03(1);
    lenpt01 = size(pt01);
    lenpt01 = lenpt01(1);
    lenpt003 = size(pt003);
    lenpt003 = lenpt003(1);

    loss = zeros(lenpt1+lenpt06+lenpt03+lenpt01+lenpt003,1);
    count = 1;
 
    %loss for strain rate pt1
    DFGRD0 = eye(3);
    DFGRD1 = eye(3);
    STATEV = zeros(18,1);
    STRESS = zeros(6,1);
    DDSDDE = zeros(6,6);
    FIRSTPIOLA = zeros(3,3);
    STATEV(1) = 1.0;
    STATEV(2) = 1.0;
    STATEV(3) = 1.0;
    STATEV(4) = 0.0;
    STATEV(5) = 0.0;
    STATEV(6) = 0.0;
    STATEV(7) = 1.0;
    STATEV(8) = 1.0;
    STATEV(9) = 1.0;
    STATEV(13) = 1.0;
    STATEV(14) = 1.0;
    STATEV(15) = 1.0;

    for datapoint = 2:lenpt1
        STRETCH_OLD = pt1(datapoint-1, 2);
        STRETCH_NEW = pt1(datapoint, 2);
        DTIME = pt1(datapoint, 4) - pt1(datapoint-1, 4);
    
        DFGRD0(1,1) = STRETCH_OLD;
        DFGRD0(2,2) = STRETCH_OLD^(-.5);
        DFGRD0(3,3) = STRETCH_OLD^(-.5);
    
        DFGRD1(1,1) = STRETCH_NEW;
        DFGRD1(2,2) = STRETCH_NEW^(-.5);
        DFGRD1(3,3) = STRETCH_NEW^(-.5);
        [STATEV, STRESS, DDSDDE, FIRSTPIOLA, exitflag] = VISCOEL_OGDEN_3EL(DFGRD0, DFGRD1, STATEV, DTIME, ...
                                                             MU, ALPHA, KELAS, ...
                                                             MUVIS, ALPHAVIS, KVIS, ETADEV, ETAVOL, ...
                                                             MUVIS_2, ALPHAVIS_2, KVIS_2, ETADEV_2, ETAVOL_2, ...
                                                             MUVIS_3, ALPHAVIS_3, KVIS_3, ETADEV_3, ETAVOL_3);
        loss(count) = (FIRSTPIOLA(1,1)*1000.0 - pt1(datapoint,1)); % Stress Data in kPa
        FIRSTPIOLA(2,2)
        count = count+1;
    end
    
 % loss for strain rate pt06
    DFGRD0 = eye(3);
    DFGRD1 = eye(3);
    STATEV = zeros(18,1);
    STRESS = zeros(6,1);
    DDSDDE = zeros(6,6);
    FIRSTPIOLA = zeros(3,3);
    STATEV(1) = 1.0;
    STATEV(2) = 1.0;
    STATEV(3) = 1.0;
    STATEV(4) = 0.0;
    STATEV(5) = 0.0;
    STATEV(6) = 0.0;
    STATEV(7) = 1.0;
    STATEV(8) = 1.0;
    STATEV(9) = 1.0;
    STATEV(13) = 1.0;
    STATEV(14) = 1.0;
    STATEV(15) = 1.0;
    
    for datapoint = 2:lenpt06
        STRETCH_OLD = pt06(datapoint-1, 2);
        STRETCH_NEW = pt06(datapoint, 2);
        DTIME = pt06(datapoint, 4) - pt06(datapoint-1, 4);
    
    
        DFGRD0(1,1) = STRETCH_OLD;
        DFGRD0(2,2) = STRETCH_OLD^(-.5);
        DFGRD0(3,3) = STRETCH_OLD^(-.5);
    
        DFGRD1(1,1) = STRETCH_NEW;
        DFGRD1(2,2) = STRETCH_NEW^(-.5);
        DFGRD1(3,3) = STRETCH_NEW^(-.5);
        [STATEV, STRESS, DDSDDE, FIRSTPIOLA, exitflag] = VISCOEL_OGDEN_3EL(DFGRD0, DFGRD1, STATEV, DTIME, ...
                                                       MU, ALPHA, KELAS, ...
                                                       MUVIS, ALPHAVIS, KVIS, ETADEV, ETAVOL, ...
                                                       MUVIS_2, ALPHAVIS_2, KVIS_2, ETADEV_2, ETAVOL_2, ...
                                                       MUVIS_3, ALPHAVIS_3, KVIS_3, ETADEV_3, ETAVOL_3);
        loss(count) = (FIRSTPIOLA(1,1)*1000.0 - pt06(datapoint,1)); % Stress Data in kPa
        count = count+1;
    end
    
    %loss for strain rate pt03
    DFGRD0 = eye(3);
    DFGRD1 = eye(3);
    STATEV = zeros(18,1);
    STRESS = zeros(6,1);
    DDSDDE = zeros(6,6);
    FIRSTPIOLA = zeros(3,3);
    STATEV(1) = 1.0;
    STATEV(2) = 1.0;
    STATEV(3) = 1.0;
    STATEV(4) = 0.0;
    STATEV(5) = 0.0;
    STATEV(6) = 0.0;
    STATEV(7) = 1.0;
    STATEV(8) = 1.0;
    STATEV(9) = 1.0;
    STATEV(13) = 1.0;
    STATEV(14) = 1.0;
    STATEV(15) = 1.0;
    
    for datapoint = 2:lenpt03
        STRETCH_OLD = pt03(datapoint-1, 2);
        STRETCH_NEW = pt03(datapoint, 2);
        DTIME = pt03(datapoint, 4) - pt03(datapoint-1, 4);
    
    
        DFGRD0(1,1) = STRETCH_OLD;
        DFGRD0(2,2) = STRETCH_OLD^(-.5);
        DFGRD0(3,3) = STRETCH_OLD^(-.5);
    
        DFGRD1(1,1) = STRETCH_NEW;
        DFGRD1(2,2) = STRETCH_NEW^(-.5);
        DFGRD1(3,3) = STRETCH_NEW^(-.5);
        [STATEV, STRESS, DDSDDE, FIRSTPIOLA, exitflag] = VISCOEL_OGDEN_3EL(DFGRD0, DFGRD1, STATEV, DTIME, ...
                                                             MU, ALPHA, KELAS, ...
                                                             MUVIS, ALPHAVIS, KVIS, ETADEV, ETAVOL, ...
                                                             MUVIS_2, ALPHAVIS_2, KVIS_2, ETADEV_2, ETAVOL_2, ...
                                                             MUVIS_3, ALPHAVIS_3, KVIS_3, ETADEV_3, ETAVOL_3);
        loss(count) = (FIRSTPIOLA(1,1)*1000.0 - pt03(datapoint,1)); % Stress Data in kPa
        count = count+1; 
    end
    
  % loss for strain rate pt01
    DFGRD0 = eye(3);
    DFGRD1 = eye(3);
    STATEV = zeros(18,1);
    STRESS = zeros(6,1);
    DDSDDE = zeros(6,6);
    FIRSTPIOLA = zeros(3,3);
    STATEV(1) = 1.0;
    STATEV(2) = 1.0;
    STATEV(3) = 1.0;
    STATEV(4) = 0.0;
    STATEV(5) = 0.0;
    STATEV(6) = 0.0;
    STATEV(7) = 1.0;
    STATEV(8) = 1.0;
    STATEV(9) = 1.0;
    STATEV(13) = 1.0;
    STATEV(14) = 1.0;
    STATEV(15) = 1.0;
    
    for datapoint = 2:lenpt01
        STRETCH_OLD = pt01(datapoint-1, 2);
        STRETCH_NEW = pt01(datapoint, 2);
        DTIME = pt01(datapoint, 4) - pt01(datapoint-1, 4);
    
    
        DFGRD0(1,1) = STRETCH_OLD;
        DFGRD0(2,2) = STRETCH_OLD^(-.5);
        DFGRD0(3,3) = STRETCH_OLD^(-.5);
    
        DFGRD1(1,1) = STRETCH_NEW;
        DFGRD1(2,2) = STRETCH_NEW^(-.5);
        DFGRD1(3,3) = STRETCH_NEW^(-.5);
        [STATEV, STRESS, DDSDDE, FIRSTPIOLA, exitflag] = VISCOEL_OGDEN_3EL(DFGRD0, DFGRD1, STATEV, DTIME, ...
                                                       MU, ALPHA, KELAS, ...
                                                       MUVIS, ALPHAVIS, KVIS, ETADEV, ETAVOL, ...
                                                       MUVIS_2, ALPHAVIS_2, KVIS_2, ETADEV_2, ETAVOL_2, ...
                                                       MUVIS_3, ALPHAVIS_3, KVIS_3, ETADEV_3, ETAVOL_3);
        loss(count) = (FIRSTPIOLA(1,1)*1000.0 - pt01(datapoint,1)); % Stress Data in kPa
        count = count+1;
    end
    
    % loss for strain rate pt003
    DFGRD0 = eye(3);
    DFGRD1 = eye(3);
    STATEV = zeros(18,1);
    STRESS = zeros(6,1);
    DDSDDE = zeros(6,6);
    FIRSTPIOLA = zeros(3,3);
    STATEV(1) = 1.0;
    STATEV(2) = 1.0;
    STATEV(3) = 1.0;
    STATEV(4) = 0.0;
    STATEV(5) = 0.0;
    STATEV(6) = 0.0;
    STATEV(7) = 1.0;
    STATEV(8) = 1.0;
    STATEV(9) = 1.0;
    STATEV(13) = 1.0;
    STATEV(14) = 1.0;
    STATEV(15) = 1.0;
    
    for datapoint = 2:lenpt003
        STRETCH_OLD = pt003(datapoint-1, 2);
        STRETCH_NEW = pt003(datapoint, 2);
        DTIME = pt003(datapoint, 4) - pt003(datapoint-1, 4);
    
    
        DFGRD0(1,1) = STRETCH_OLD;
        DFGRD0(2,2) = STRETCH_OLD^(-.5);
        DFGRD0(3,3) = STRETCH_OLD^(-.5);
    
        DFGRD1(1,1) = STRETCH_NEW;
        DFGRD1(2,2) = STRETCH_NEW^(-.5);
        DFGRD1(3,3) = STRETCH_NEW^(-.5);
        [STATEV, STRESS, DDSDDE, FIRSTPIOLA, exitflag] = VISCOEL_OGDEN_3EL(DFGRD0, DFGRD1, STATEV, DTIME, ...
                                                       MU, ALPHA, KELAS, ...
                                                       MUVIS, ALPHAVIS, KVIS, ETADEV, ETAVOL, ...
                                                       MUVIS_2, ALPHAVIS_2, KVIS_2, ETADEV_2, ETAVOL_2, ...
                                                       MUVIS_3, ALPHAVIS_3, KVIS_3, ETADEV_3, ETAVOL_3);
        loss(count) = (FIRSTPIOLA(1,1)*1000.0 - pt003(datapoint,1)); % Stress Data in kPa
        count = count+1;
    end
end