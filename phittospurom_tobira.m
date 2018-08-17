function [ala,hot,refl,tran]=phittospurom_tobira(lambda)

% Mar�a Eugenia Beget1 and V�ctor Bettachini2 
% 1 Los Reseros y Las Caba�as s/n, 1686 Hurlingham, Buenos Aires, Argentina. Tel.: +54 11 4621 1684; fax: +54 11 4621 5663.
% mbeget@cnia.inta.gov.ar 
% 2 Laboratorio de Optoelectr�nica, Instituto Tecnol�gico de Buenos Aires, Av. Eduardo Madero 399, C1106ACD Buenos Aires, Argentina.
% vbettach@itba.edu.ar

% Beget, M.E., V. A. Bettachini, C. Di Bella, F. Baret. 2013.  
% SAILHFlood: a radiative transfer model for flooded vegetation. Ecological Modelling, 257:25-35.




% Phittospurom tobira

% Leaf optical properties meassured at Unit� Climat - INRA Avignon -
% France, March 2007
% M.E. Beget

format short e;

% INPUTS: ala, hot, Tablerefl, Tabletran

ala=39;         % ala  : mean leaf inclination angle (deegres)
hot=0.42;       % hot : hot spot parameter
% Tablerefl : leaf reflectance at 1 nm resolution (400 - 2400 nm)
Tablerefl=[0.0468;0.0477;0.0485;0.0493;0.0501;0.0508;0.0515;0.0521;0.0527;0.0533;0.0538;0.0543;0.0549;0.0554;0.0559;0.0564;0.0569;0.0574;0.0579;0.0584;0.0588;0.0592;0.0596;0.06;0.0603;0.0606;0.0608;0.061;0.0612;0.0614;0.0616;0.0617;0.0618;0.0619;0.062;0.0622;0.0623;0.0624;0.0626;0.0629;0.0631;0.0634;0.0638;0.0641;0.0645;0.0649;0.0654;0.0658;0.0662;0.0667;0.0671;0.0676;0.068;0.0684;0.0688;0.0691;0.0695;0.0697;0.07;0.0702;0.0704;0.0705;0.0707;0.0708;0.0708;0.0709;0.0709;0.071;0.071;0.071;0.0711;0.0711;0.0711;0.0712;0.0712;0.0713;0.0714;0.0715;0.0715;0.0716;0.0717;0.0719;0.072;0.0721;0.0723;0.0724;0.0726;0.0729;0.0731;0.0734;0.0738;0.0742;0.0747;0.0752;0.0758;0.0765;0.0772;0.078;0.0789;0.0799;0.081;0.0821;0.0833;0.0845;0.0859;0.0873;0.0888;0.0904;0.092;0.0937;0.0955;0.0974;0.0994;0.1014;0.1035;0.1057;0.1079;0.1101;0.1124;0.1147;0.117;0.1193;0.1216;0.1238;0.1261;0.1283;0.1304;0.1324;0.1344;0.1362;0.138;0.1396;0.1411;0.1426;0.1439;0.1451;0.1462;0.1472;0.1481;0.1489;0.1497;0.1505;0.1511;0.1518;0.1524;0.1529;0.1534;0.1539;0.1544;0.1548;0.1551;0.1553;0.1555;0.1556;0.1556;0.1555;0.1552;0.1548;0.1544;0.1537;0.153;0.1522;0.1513;0.1503;0.1492;0.1481;0.1469;0.1457;0.1445;0.1432;0.142;0.1407;0.1394;0.1382;0.137;0.1358;0.1346;0.1335;0.1325;0.1315;0.1305;0.1296;0.1288;0.128;0.1273;0.1266;0.1259;0.1253;0.1247;0.1242;0.1237;0.1232;0.1228;0.1224;0.122;0.1217;0.1213;0.121;0.1207;0.1204;0.1201;0.1198;0.1195;0.1191;0.1187;0.1182;0.1177;0.1172;0.1166;0.116;0.1154;0.1147;0.114;0.1134;0.1127;0.112;0.1114;0.1108;0.1102;0.1097;0.1092;0.1087;0.1084;0.1081;0.1078;0.1076;0.1074;0.1072;0.107;0.1069;0.1067;0.1065;0.1062;0.1059;0.1055;0.105;0.1044;0.1037;0.1029;0.1021;0.1011;0.1001;0.0991;0.098;0.097;0.0959;0.0949;0.0939;0.093;0.0922;0.0914;0.0906;0.0899;0.0892;0.0886;0.0879;0.0873;0.0866;0.0859;0.0851;0.0843;0.0835;0.0827;0.0819;0.081;0.0802;0.0794;0.0787;0.078;0.0774;0.0769;0.0765;0.0761;0.0759;0.0758;0.0759;0.076;0.0764;0.0769;0.0776;0.0786;0.0798;0.0812;0.083;0.0851;0.0875;0.0903;0.0935;0.0971;0.1011;0.1055;0.1103;0.1155;0.1211;0.127;0.1333;0.1399;0.1469;0.1541;0.1615;0.1692;0.177;0.185;0.193;0.2012;0.2094;0.2176;0.2259;0.2342;0.2424;0.2507;0.259;0.2672;0.2754;0.2836;0.2918;0.2999;0.308;0.3159;0.3238;0.3316;0.3392;0.3467;0.354;0.3611;0.3681;0.3747;0.3812;0.3874;0.3933;0.399;0.4044;0.4095;0.4143;0.4189;0.4232;0.4272;0.431;0.4344;0.4377;0.4407;0.4434;0.4459;0.4482;0.4503;0.4522;0.4539;0.4555;0.4569;0.4581;0.4592;0.4602;0.4611;0.4618;0.4625;0.4631;0.4636;0.4641;0.4645;0.4648;0.4651;0.4653;0.4655;0.4657;0.4658;0.4659;0.466;0.4661;0.4661;0.4662;0.4662;0.4662;0.4662;0.4662;0.4662;0.4662;0.4663;0.4663;0.4663;0.4663;0.4663;0.4664;0.4664;0.4664;0.4665;0.4665;0.4666;0.4666;0.4667;0.4667;0.4668;0.4668;0.4669;0.467;0.467;0.4671;0.4672;0.4673;0.4673;0.4674;0.4675;0.4676;0.4676;0.4677;0.4678;0.4678;0.4679;0.4679;0.468;0.4681;0.4681;0.4682;0.4682;0.4683;0.4684;0.4685;0.4685;0.4686;0.4687;0.4688;0.4688;0.4689;0.469;0.469;0.4691;0.4692;0.4692;0.4693;0.4694;0.4694;0.4695;0.4695;0.4696;0.4696;0.4697;0.4697;0.4698;0.4699;0.4699;0.47;0.47;0.4701;0.4701;0.4702;0.4702;0.4702;0.4703;0.4703;0.4704;0.4704;0.4705;0.4706;0.4706;0.4707;0.4707;0.4708;0.4708;0.4709;0.4709;0.471;0.471;0.4711;0.4711;0.4711;0.4711;0.4712;0.4712;0.4712;0.4712;0.4712;0.4711;0.4711;0.4711;0.4711;0.4711;0.4712;0.4712;0.4712;0.4712;0.4712;0.4712;0.4712;0.4712;0.4712;0.4711;0.4711;0.4711;0.4711;0.471;0.471;0.471;0.4709;0.4709;0.4708;0.4708;0.4707;0.4707;0.4706;0.4705;0.4705;0.4705;0.4704;0.4704;0.4703;0.4703;0.4702;0.4702;0.4701;0.4701;0.47;0.47;0.4699;0.4698;0.4698;0.4697;0.4696;0.4696;0.4695;0.4695;0.4695;0.4694;0.4694;0.4693;0.4692;0.4691;0.469;0.4689;0.4687;0.4686;0.4684;0.4682;0.468;0.4678;0.4676;0.4674;0.4672;0.467;0.4668;0.4665;0.4663;0.4661;0.4658;0.4656;0.4653;0.465;0.4646;0.4643;0.4639;0.4635;0.4631;0.4627;0.4623;0.4618;0.4614;0.4609;0.4604;0.46;0.4595;0.4591;0.4586;0.4583;0.4579;0.4576;0.4573;0.4571;0.4569;0.4567;0.4565;0.4564;0.4562;0.4561;0.4559;0.4558;0.4556;0.4555;0.4554;0.4554;0.4555;0.4556;0.4557;0.4558;0.4559;0.456;0.4562;0.4563;0.4564;0.4565;0.4566;0.4568;0.4569;0.457;0.4571;0.4573;0.4574;0.4575;0.4577;0.4578;0.458;0.4582;0.4583;0.4585;0.4587;0.4589;0.4591;0.4593;0.4596;0.4598;0.46;0.4603;0.4605;0.4607;0.461;0.4612;0.4614;0.4617;0.4619;0.4622;0.4624;0.4626;0.4629;0.4631;0.4633;0.4636;0.4638;0.464;0.4642;0.4644;0.4646;0.4648;0.4649;0.4651;0.4653;0.4654;0.4656;0.4657;0.4658;0.466;0.4661;0.4662;0.4663;0.4664;0.4665;0.4666;0.4667;0.4668;0.4669;0.4669;0.467;0.4671;0.4671;0.4672;0.4672;0.4672;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4673;0.4672;0.4672;0.4672;0.4672;0.4672;0.4671;0.4671;0.4671;0.4671;0.467;0.467;0.4669;0.4669;0.4668;0.4668;0.4667;0.4667;0.4666;0.4665;0.4665;0.4664;0.4663;0.4662;0.4661;0.4661;0.466;0.4659;0.4658;0.4656;0.4655;0.4654;0.4652;0.465;0.4648;0.4645;0.4642;0.4639;0.4636;0.4632;0.4627;0.4622;0.4617;0.4611;0.4604;0.4598;0.459;0.4582;0.4574;0.4565;0.4556;0.4546;0.4536;0.4525;0.4514;0.4503;0.4491;0.448;0.4468;0.4456;0.4444;0.4432;0.4419;0.4408;0.4396;0.4384;0.4373;0.4362;0.4351;0.4341;0.4331;0.4321;0.4312;0.4303;0.4295;0.4287;0.428;0.4273;0.4267;0.4261;0.4255;0.425;0.4245;0.424;0.4236;0.4232;0.4229;0.4225;0.4222;0.4219;0.4216;0.4213;0.421;0.4208;0.4205;0.4202;0.42;0.4197;0.4194;0.4191;0.4189;0.4186;0.4183;0.4181;0.4178;0.4175;0.4173;0.4171;0.4169;0.4167;0.4165;0.4164;0.4162;0.4162;0.4161;0.416;0.416;0.416;0.4161;0.4162;0.4162;0.4164;0.4165;0.4167;0.4169;0.4171;0.4173;0.4176;0.4179;0.4181;0.4185;0.4188;0.4191;0.4195;0.4198;0.4202;0.4206;0.421;0.4214;0.4218;0.4222;0.4226;0.4231;0.4235;0.4239;0.4243;0.4247;0.4251;0.4255;0.4259;0.4262;0.4266;0.4269;0.4272;0.4275;0.4278;0.4281;0.4283;0.4286;0.4288;0.429;0.4292;0.4294;0.4295;0.4297;0.4298;0.43;0.4301;0.4302;0.4303;0.4304;0.4305;0.4306;0.4307;0.4308;0.4308;0.4309;0.4309;0.431;0.431;0.431;0.431;0.431;0.4309;0.4309;0.4309;0.4308;0.4307;0.4306;0.4305;0.4304;0.4303;0.4302;0.43;0.4298;0.4297;0.4295;0.4293;0.4291;0.4289;0.4286;0.4284;0.4281;0.4279;0.4276;0.4273;0.427;0.4267;0.4264;0.4261;0.4257;0.4254;0.425;0.4246;0.4242;0.4237;0.4233;0.4228;0.4223;0.4217;0.4211;0.4205;0.4199;0.4192;0.4185;0.4178;0.417;0.4163;0.4155;0.4146;0.4138;0.4129;0.4119;0.411;0.41;0.409;0.408;0.407;0.4059;0.4048;0.4037;0.4025;0.4014;0.4002;0.3989;0.3977;0.3965;0.3952;0.394;0.3927;0.3914;0.3902;0.3889;0.3877;0.3865;0.3853;0.3841;0.3829;0.3817;0.3806;0.3794;0.3783;0.3772;0.3761;0.375;0.374;0.3729;0.3718;0.3708;0.3697;0.3687;0.3676;0.3665;0.3654;0.3642;0.363;0.3618;0.3605;0.3592;0.3578;0.3563;0.3547;0.353;0.3512;0.3493;0.3473;0.3452;0.3429;0.3405;0.3379;0.3352;0.3323;0.3293;0.3261;0.3227;0.3191;0.3154;0.3115;0.3075;0.3033;0.2989;0.2945;0.2899;0.2852;0.2804;0.2756;0.2707;0.2658;0.2608;0.2559;0.251;0.2461;0.2413;0.2365;0.2318;0.2272;0.2228;0.2185;0.2142;0.2102;0.2063;0.2025;0.1989;0.1955;0.1922;0.1891;0.1861;0.1833;0.1806;0.1781;0.1757;0.1735;0.1714;0.1695;0.1676;0.1659;0.1643;0.1628;0.1614;0.1601;0.1589;0.1578;0.1568;0.1558;0.155;0.1542;0.1534;0.1528;0.1521;0.1516;0.1511;0.1506;0.1502;0.1498;0.1495;0.1491;0.1489;0.1486;0.1484;0.1483;0.1481;0.148;0.148;0.1479;0.1479;0.1479;0.1479;0.148;0.1481;0.1481;0.1483;0.1484;0.1486;0.1487;0.1489;0.1491;0.1494;0.1496;0.1499;0.1503;0.1506;0.151;0.1514;0.1519;0.1524;0.1529;0.1535;0.1541;0.1548;0.1554;0.1561;0.1569;0.1576;0.1584;0.1592;0.16;0.1608;0.1616;0.1624;0.1633;0.1642;0.165;0.1659;0.1668;0.1678;0.1687;0.1697;0.1706;0.1716;0.1726;0.1737;0.1747;0.1757;0.1768;0.1779;0.1789;0.18;0.1811;0.1821;0.1832;0.1843;0.1853;0.1864;0.1874;0.1885;0.1896;0.1906;0.1917;0.1928;0.1938;0.1949;0.196;0.1971;0.1982;0.1993;0.2004;0.2015;0.2026;0.2037;0.2048;0.2058;0.2069;0.2079;0.2089;0.2099;0.2109;0.2119;0.2128;0.2138;0.2147;0.2157;0.2166;0.2176;0.2186;0.2195;0.2205;0.2215;0.2225;0.2234;0.2244;0.2254;0.2263;0.2273;0.2282;0.2291;0.23;0.2309;0.2318;0.2326;0.2334;0.2342;0.235;0.2358;0.2366;0.2374;0.2382;0.2389;0.2397;0.2405;0.2412;0.242;0.2428;0.2436;0.2443;0.2451;0.2458;0.2466;0.2473;0.248;0.2487;0.2494;0.2501;0.2508;0.2515;0.2521;0.2528;0.2535;0.2541;0.2548;0.2554;0.2561;0.2568;0.2574;0.2581;0.2587;0.2594;0.26;0.2607;0.2613;0.2619;0.2625;0.2631;0.2637;0.2642;0.2648;0.2654;0.2659;0.2665;0.267;0.2676;0.2681;0.2687;0.2692;0.2698;0.2704;0.2709;0.2715;0.2721;0.2726;0.2732;0.2737;0.2743;0.2748;0.2753;0.2759;0.2764;0.2769;0.2774;0.2779;0.2784;0.2789;0.2794;0.2799;0.2803;0.2808;0.2813;0.2817;0.2822;0.2826;0.2831;0.2835;0.2839;0.2843;0.2847;0.2851;0.2855;0.2858;0.2862;0.2865;0.2869;0.2872;0.2875;0.2879;0.2882;0.2885;0.2888;0.2891;0.2894;0.2897;0.29;0.2902;0.2904;0.2906;0.2908;0.291;0.2911;0.2912;0.2913;0.2914;0.2914;0.2915;0.2915;0.2915;0.2914;0.2914;0.2913;0.2912;0.2911;0.2909;0.2907;0.2905;0.2903;0.29;0.2896;0.2893;0.2889;0.2884;0.288;0.2875;0.2869;0.2864;0.2858;0.2852;0.2846;0.2839;0.2833;0.2827;0.282;0.2814;0.2807;0.28;0.2794;0.2787;0.278;0.2774;0.2767;0.276;0.2753;0.2746;0.2739;0.2733;0.2726;0.2719;0.2713;0.2707;0.27;0.2694;0.2688;0.2682;0.2676;0.267;0.2664;0.2659;0.2653;0.2648;0.2642;0.2637;0.2631;0.2626;0.2622;0.2617;0.2613;0.2609;0.2605;0.2602;0.2599;0.2596;0.2593;0.2591;0.2589;0.2588;0.2586;0.2584;0.2582;0.2581;0.2579;0.2577;0.2575;0.2572;0.257;0.2567;0.2564;0.2561;0.2558;0.2554;0.2551;0.2547;0.2543;0.254;0.2536;0.2531;0.2527;0.2523;0.2519;0.2514;0.251;0.2505;0.2501;0.2496;0.2492;0.2488;0.2483;0.2479;0.2475;0.2472;0.2468;0.2464;0.2461;0.2457;0.2454;0.2452;0.245;0.2449;0.2449;0.2449;0.245;0.2451;0.2452;0.2453;0.2454;0.2455;0.2457;0.2458;0.2459;0.246;0.2461;0.2463;0.2464;0.2465;0.2467;0.2468;0.247;0.2471;0.2473;0.2474;0.2476;0.2478;0.2479;0.2481;0.2483;0.2484;0.2486;0.2488;0.249;0.2492;0.2493;0.2495;0.2497;0.2499;0.2501;0.2502;0.2504;0.2506;0.2507;0.2509;0.251;0.2511;0.2512;0.2513;0.2514;0.2514;0.2514;0.2514;0.2513;0.2512;0.2511;0.2509;0.2507;0.2504;0.2501;0.2498;0.2493;0.2489;0.2483;0.2477;0.2471;0.2464;0.2455;0.2447;0.2437;0.2427;0.2415;0.2403;0.239;0.2376;0.2361;0.2345;0.2329;0.2311;0.2292;0.2271;0.225;0.2228;0.2205;0.218;0.2155;0.2128;0.21;0.2072;0.2042;0.2011;0.198;0.1947;0.1914;0.1879;0.1844;0.1809;0.1772;0.1736;0.1698;0.1661;0.1623;0.1584;0.1546;0.1508;0.1469;0.1431;0.1393;0.1355;0.1318;0.1281;0.1245;0.1209;0.1174;0.114;0.1107;0.1074;0.1043;0.1013;0.0983;0.0955;0.0928;0.0902;0.0878;0.0854;0.0832;0.0811;0.0791;0.0772;0.0755;0.0738;0.0723;0.0709;0.0695;0.0683;0.0672;0.0662;0.0652;0.0644;0.0637;0.063;0.0624;0.0619;0.0615;0.0611;0.0608;0.0605;0.0603;0.0602;0.0601;0.06;0.06;0.06;0.0601;0.0601;0.0603;0.0604;0.0605;0.0607;0.0609;0.0611;0.0613;0.0615;0.0618;0.0621;0.0623;0.0626;0.0629;0.0631;0.0634;0.0637;0.064;0.0643;0.0645;0.0648;0.0651;0.0654;0.0657;0.0659;0.0662;0.0665;0.0668;0.0671;0.0674;0.0677;0.068;0.0683;0.0686;0.0689;0.0693;0.0696;0.07;0.0703;0.0707;0.071;0.0714;0.0718;0.0721;0.0725;0.0729;0.0733;0.0737;0.0741;0.0745;0.0749;0.0753;0.0757;0.0761;0.0766;0.077;0.0775;0.078;0.0785;0.079;0.0795;0.08;0.0805;0.0811;0.0816;0.0821;0.0827;0.0832;0.0838;0.0843;0.0849;0.0855;0.086;0.0866;0.0871;0.0877;0.0883;0.0888;0.0893;0.0899;0.0904;0.091;0.0915;0.092;0.0925;0.093;0.0935;0.0941;0.0945;0.095;0.0955;0.096;0.0965;0.097;0.0974;0.0979;0.0984;0.0989;0.0993;0.0998;0.1003;0.1007;0.1012;0.1016;0.1021;0.1025;0.103;0.1035;0.1039;0.1044;0.1048;0.1053;0.1057;0.1062;0.1066;0.1071;0.1076;0.108;0.1085;0.1089;0.1094;0.1098;0.1103;0.1107;0.1112;0.1116;0.1121;0.1125;0.113;0.1134;0.1138;0.1143;0.1147;0.1151;0.1155;0.116;0.1164;0.1168;0.1172;0.1177;0.1181;0.1185;0.1189;0.1194;0.1198;0.1202;0.1206;0.121;0.1214;0.1218;0.1222;0.1226;0.1229;0.1233;0.1237;0.124;0.1244;0.1247;0.1251;0.1254;0.1257;0.1261;0.1264;0.1267;0.1271;0.1274;0.1277;0.128;0.1283;0.1286;0.1289;0.1291;0.1294;0.1296;0.1299;0.1301;0.1304;0.1306;0.1308;0.1311;0.1313;0.1315;0.1318;0.132;0.1323;0.1326;0.1329;0.1332;0.1335;0.1339;0.1342;0.1346;0.135;0.1354;0.1358;0.1362;0.1366;0.1371;0.1375;0.1379;0.1383;0.1387;0.1392;0.1396;0.14;0.1404;0.1409;0.1413;0.1417;0.1421;0.1426;0.143;0.1434;0.1439;0.1443;0.1448;0.1452;0.1457;0.1461;0.1466;0.1471;0.1475;0.148;0.1484;0.1489;0.1493;0.1498;0.1502;0.1506;0.1511;0.1515;0.1519;0.1523;0.1527;0.1532;0.1536;0.154;0.1544;0.1548;0.1552;0.1556;0.156;0.1564;0.1568;0.1571;0.1575;0.1579;0.1582;0.1586;0.1589;0.1592;0.1595;0.1598;0.1601;0.1604;0.1607;0.1609;0.1612;0.1614;0.1616;0.1618;0.162;0.1622;0.1624;0.1625;0.1627;0.1628;0.1629;0.163;0.1631;0.1631;0.1632;0.1633;0.1633;0.1633;0.1634;0.1634;0.1634;0.1634;0.1634;0.1634;0.1634;0.1635;0.1634;0.1634;0.1634;0.1634;0.1633;0.1633;0.1632;0.163;0.1629;0.1627;0.1626;0.1624;0.1621;0.1619;0.1617;0.1614;0.1611;0.1608;0.1606;0.1603;0.16;0.1597;0.1594;0.1592;0.1589;0.1586;0.1583;0.1579;0.1576;0.1572;0.1568;0.1564;0.156;0.1555;0.155;0.1545;0.1539;0.1533;0.1528;0.1521;0.1515;0.1509;0.1502;0.1496;0.149;0.1483;0.1477;0.1471;0.1465;0.1459;0.1453;0.1447;0.1441;0.1435;0.143;0.1424;0.1418;0.1412;0.1407;0.1401;0.1395;0.139;0.1384;0.1379;0.1373;0.1367;0.1362;0.1356;0.1351;0.1345;0.134;0.1334;0.1329;0.1323;0.1318;0.1313;0.1307;0.1302;0.1296;0.1291;0.1285;0.128;0.1275;0.1269;0.1264;0.1259;0.1254;0.1249;0.1244;0.1239;0.1235;0.1231;0.1226;0.1222;0.1219;0.1215;0.1212;0.1209;0.1205;0.1203;0.12;0.1197;0.1195;0.1192;0.119;0.1188;0.1185;0.1183;0.1181;0.1179;0.1177;0.1175;0.1174;0.1172;0.117;0.1168;0.1167;0.1165;0.1163;0.1161;0.1158;0.1156;0.1153;0.115;0.1147;0.1144;0.114;0.1137;0.1133;0.1129;0.1125;0.1121;0.1117;0.1113;0.1108;0.1104;0.11;0.1096;0.1091;0.1087;0.1082;0.1077;0.1072;0.1067;0.1061;0.1055;0.1049;0.1042;0.1036;0.1029;0.1022;0.1016;0.1009;0.1003;0.0997;0.099;0.0985;0.0979;0.0974;0.0969;0.0964;0.096;0.0956;0.0952;0.0949;0.0946;0.0944;0.0942;0.094;0.0939;0.0938;0.0939;0.094;0.0941;0.0944;0.0947;0.0951;0.0956;0.0962;0.0969;0.0976;0.0984;0.0992;0.1001;0.101;0.1019;0.1028;0.1037;0.1047;0.1056];           % Tablerefl, 2001 coefficients (400 nm to 2400 nm)  
% Tabletran : leaf transmittance at 1 nm resolution (400 - 2400 nm)
Tabletran=[0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0002;0.0002;0.0002;0.0002;0.0002;0.0002;0.0002;0.0002;0.0002;0.0003;0.0003;0.0003;0.0003;0.0003;0.0004;0.0004;0.0004;0.0004;0.0004;0.0005;0.0005;0.0005;0.0005;0.0005;0.0006;0.0006;0.0007;0.0007;0.0008;0.0008;0.0009;0.001;0.001;0.0011;0.0011;0.0011;0.0011;0.0012;0.0012;0.0012;0.0012;0.0012;0.0012;0.0012;0.0013;0.0013;0.0014;0.0014;0.0015;0.0015;0.0016;0.0016;0.0017;0.0017;0.0017;0.0017;0.0017;0.0017;0.0017;0.0018;0.0018;0.0018;0.0018;0.0019;0.0019;0.002;0.002;0.0021;0.0022;0.0022;0.0023;0.0024;0.0025;0.0026;0.0027;0.0029;0.003;0.0032;0.0033;0.0035;0.0038;0.004;0.0043;0.0046;0.005;0.0054;0.0058;0.0063;0.0069;0.0075;0.0081;0.0089;0.0097;0.0105;0.0115;0.0124;0.0135;0.0145;0.0157;0.0168;0.0179;0.0191;0.0203;0.0215;0.0226;0.0238;0.0249;0.0259;0.027;0.028;0.0289;0.0298;0.0306;0.0314;0.0321;0.0327;0.0334;0.0339;0.0345;0.035;0.0355;0.0359;0.0364;0.0368;0.0372;0.0375;0.0379;0.0382;0.0384;0.0386;0.0387;0.0387;0.0386;0.0384;0.0381;0.0378;0.0373;0.0368;0.0362;0.0356;0.0349;0.0342;0.0335;0.0327;0.0319;0.0311;0.0303;0.0295;0.0287;0.028;0.0272;0.0265;0.0258;0.0251;0.0245;0.0239;0.0234;0.0229;0.0224;0.022;0.0216;0.0212;0.0209;0.0206;0.0203;0.02;0.0198;0.0195;0.0193;0.0192;0.019;0.0189;0.0188;0.0187;0.0186;0.0184;0.0183;0.0182;0.0181;0.0179;0.0177;0.0175;0.0173;0.017;0.0168;0.0165;0.0162;0.0159;0.0156;0.0153;0.015;0.0147;0.0145;0.0142;0.014;0.0138;0.0136;0.0134;0.0132;0.0131;0.013;0.0129;0.0129;0.0128;0.0128;0.0127;0.0127;0.0126;0.0126;0.0125;0.0123;0.0122;0.012;0.0118;0.0115;0.0112;0.0109;0.0106;0.0103;0.01;0.0097;0.0094;0.0091;0.0088;0.0086;0.0083;0.0081;0.0079;0.0077;0.0076;0.0074;0.0073;0.0071;0.0069;0.0068;0.0066;0.0065;0.0063;0.0061;0.006;0.0058;0.0057;0.0055;0.0054;0.0052;0.0051;0.0049;0.0048;0.0047;0.0046;0.0045;0.0044;0.0044;0.0044;0.0043;0.0044;0.0044;0.0045;0.0046;0.0049;0.0052;0.0057;0.0064;0.0073;0.0084;0.0099;0.0118;0.014;0.0167;0.0199;0.0236;0.0278;0.0326;0.0378;0.0435;0.0497;0.0562;0.0632;0.0704;0.0779;0.0856;0.0935;0.1016;0.1097;0.1179;0.1262;0.1345;0.1428;0.1511;0.1594;0.1678;0.176;0.1843;0.1925;0.2006;0.2086;0.2166;0.2244;0.2321;0.2396;0.2469;0.254;0.2609;0.2676;0.274;0.2802;0.286;0.2917;0.297;0.302;0.3068;0.3113;0.3156;0.3196;0.3233;0.3268;0.33;0.3331;0.3359;0.3385;0.3408;0.343;0.3451;0.3469;0.3486;0.3501;0.3515;0.3528;0.3539;0.3549;0.3559;0.3568;0.3575;0.3582;0.3589;0.3595;0.36;0.3605;0.3609;0.3613;0.3617;0.362;0.3623;0.3626;0.3628;0.3631;0.3633;0.3635;0.3637;0.3639;0.3641;0.3643;0.3645;0.3647;0.3649;0.365;0.3652;0.3654;0.3655;0.3657;0.3659;0.366;0.3662;0.3664;0.3666;0.3667;0.3669;0.3671;0.3673;0.3675;0.3677;0.3679;0.3681;0.3683;0.3685;0.3687;0.3689;0.3691;0.3693;0.3694;0.3696;0.3698;0.37;0.3703;0.3705;0.3707;0.3709;0.3711;0.3713;0.3715;0.3717;0.372;0.3722;0.3724;0.3726;0.3728;0.373;0.3732;0.3734;0.3736;0.3738;0.374;0.3742;0.3744;0.3746;0.3748;0.375;0.3752;0.3754;0.3756;0.3758;0.376;0.3762;0.3763;0.3765;0.3767;0.3769;0.3771;0.3773;0.3775;0.3777;0.3778;0.378;0.3782;0.3783;0.3785;0.3786;0.3788;0.379;0.3791;0.3793;0.3795;0.3796;0.3797;0.3799;0.38;0.3801;0.3802;0.3803;0.3803;0.3804;0.3805;0.3806;0.3807;0.3807;0.3808;0.3809;0.381;0.381;0.3811;0.3812;0.3812;0.3813;0.3814;0.3815;0.3816;0.3817;0.3818;0.3819;0.3819;0.382;0.3821;0.3821;0.3821;0.3822;0.3822;0.3822;0.3822;0.3822;0.3822;0.3822;0.3822;0.3823;0.3823;0.3824;0.3824;0.3825;0.3826;0.3826;0.3827;0.3827;0.3828;0.3828;0.3828;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3829;0.3828;0.3828;0.3827;0.3826;0.3825;0.3823;0.3822;0.3821;0.382;0.3818;0.3817;0.3816;0.3814;0.3813;0.3812;0.3811;0.3809;0.3808;0.3806;0.3805;0.3803;0.3801;0.3799;0.3797;0.3795;0.3793;0.379;0.3787;0.3784;0.378;0.3776;0.3772;0.3767;0.3762;0.3757;0.3752;0.3747;0.3743;0.3738;0.3735;0.3731;0.3728;0.3726;0.3723;0.3722;0.3721;0.372;0.3719;0.3718;0.3718;0.3717;0.3717;0.3716;0.3716;0.3716;0.3716;0.3717;0.3718;0.3719;0.3719;0.372;0.3721;0.3722;0.3723;0.3724;0.3725;0.3727;0.3728;0.373;0.3732;0.3734;0.3736;0.3738;0.374;0.3742;0.3745;0.3747;0.3749;0.3752;0.3754;0.3757;0.3759;0.3762;0.3765;0.3767;0.377;0.3773;0.3776;0.3778;0.3781;0.3784;0.3787;0.379;0.3792;0.3795;0.3798;0.3801;0.3804;0.3807;0.3809;0.3812;0.3815;0.3817;0.382;0.3822;0.3824;0.3827;0.3829;0.383;0.3832;0.3834;0.3836;0.3837;0.3838;0.384;0.3841;0.3842;0.3843;0.3844;0.3845;0.3846;0.3847;0.3848;0.3849;0.385;0.385;0.3851;0.3852;0.3853;0.3854;0.3855;0.3856;0.3857;0.3857;0.3858;0.3859;0.386;0.3861;0.3862;0.3863;0.3864;0.3865;0.3866;0.3867;0.3868;0.3868;0.3869;0.387;0.3871;0.3872;0.3873;0.3873;0.3874;0.3874;0.3875;0.3875;0.3875;0.3875;0.3875;0.3875;0.3875;0.3874;0.3874;0.3874;0.3873;0.3872;0.3872;0.3871;0.387;0.387;0.3869;0.3868;0.3868;0.3867;0.3866;0.3866;0.3865;0.3865;0.3864;0.3864;0.3863;0.3863;0.3862;0.3862;0.3861;0.3861;0.386;0.3859;0.3858;0.3857;0.3856;0.3855;0.3854;0.3852;0.385;0.3848;0.3846;0.3844;0.3841;0.3838;0.3834;0.3831;0.3826;0.3822;0.3817;0.3811;0.3805;0.3799;0.3791;0.3784;0.3776;0.3767;0.3757;0.3748;0.3737;0.3726;0.3715;0.3703;0.3691;0.3679;0.3667;0.3654;0.3641;0.3629;0.3616;0.3603;0.3591;0.3579;0.3567;0.3555;0.3544;0.3533;0.3522;0.3512;0.3503;0.3494;0.3485;0.3477;0.347;0.3462;0.3456;0.3449;0.3443;0.3437;0.3432;0.3427;0.3422;0.3418;0.3413;0.3409;0.3405;0.3402;0.3398;0.3395;0.3392;0.3388;0.3385;0.3382;0.3379;0.3376;0.3373;0.337;0.3367;0.3364;0.3361;0.3358;0.3355;0.3353;0.335;0.3347;0.3345;0.3343;0.3341;0.3339;0.3337;0.3336;0.3335;0.3334;0.3334;0.3334;0.3334;0.3335;0.3336;0.3337;0.3338;0.334;0.3342;0.3345;0.3348;0.3351;0.3354;0.3358;0.3361;0.3365;0.3369;0.3374;0.3378;0.3382;0.3386;0.3391;0.3395;0.3399;0.3404;0.3408;0.3412;0.3417;0.3421;0.3425;0.3429;0.3434;0.3438;0.3442;0.3446;0.3449;0.3453;0.3456;0.346;0.3463;0.3465;0.3468;0.3471;0.3473;0.3475;0.3477;0.3479;0.3481;0.3482;0.3484;0.3485;0.3486;0.3488;0.3489;0.349;0.3491;0.3491;0.3492;0.3493;0.3494;0.3494;0.3495;0.3496;0.3496;0.3497;0.3497;0.3497;0.3497;0.3497;0.3498;0.3498;0.3497;0.3497;0.3497;0.3497;0.3496;0.3496;0.3495;0.3494;0.3493;0.3492;0.3491;0.3489;0.3487;0.3485;0.3483;0.3481;0.3479;0.3476;0.3474;0.3471;0.3468;0.3465;0.3462;0.3458;0.3455;0.3451;0.3447;0.3443;0.3438;0.3434;0.3429;0.3424;0.3419;0.3413;0.3407;0.3401;0.3395;0.3388;0.3381;0.3374;0.3366;0.3358;0.335;0.3342;0.3333;0.3324;0.3315;0.3305;0.3296;0.3286;0.3275;0.3265;0.3254;0.3243;0.3231;0.322;0.3208;0.3196;0.3183;0.3171;0.3159;0.3146;0.3133;0.3121;0.3108;0.3095;0.3083;0.307;0.3058;0.3045;0.3032;0.302;0.3008;0.2995;0.2983;0.2971;0.2959;0.2947;0.2935;0.2923;0.2911;0.29;0.2888;0.2877;0.2865;0.2854;0.2842;0.283;0.2818;0.2805;0.2792;0.2778;0.2764;0.2749;0.2733;0.2717;0.2699;0.268;0.266;0.2639;0.2617;0.2593;0.2568;0.2541;0.2512;0.2482;0.245;0.2417;0.2382;0.2345;0.2306;0.2266;0.2224;0.2181;0.2136;0.209;0.2043;0.1996;0.1947;0.1897;0.1847;0.1797;0.1747;0.1697;0.1646;0.1597;0.1548;0.1499;0.1452;0.1406;0.1361;0.1317;0.1275;0.1234;0.1195;0.1158;0.1122;0.1088;0.1056;0.1025;0.0996;0.0969;0.0944;0.092;0.0897;0.0876;0.0856;0.0838;0.0821;0.0805;0.079;0.0777;0.0764;0.0753;0.0742;0.0733;0.0724;0.0716;0.0708;0.0702;0.0696;0.069;0.0685;0.068;0.0676;0.0672;0.0669;0.0665;0.0663;0.066;0.0658;0.0656;0.0654;0.0653;0.0652;0.0651;0.065;0.065;0.065;0.065;0.0651;0.0651;0.0652;0.0653;0.0654;0.0656;0.0658;0.0659;0.0661;0.0664;0.0666;0.0669;0.0672;0.0675;0.0678;0.0682;0.0686;0.0691;0.0696;0.0701;0.0707;0.0712;0.0719;0.0725;0.0732;0.0739;0.0747;0.0754;0.0762;0.077;0.0778;0.0786;0.0795;0.0803;0.0812;0.0821;0.083;0.0839;0.0848;0.0857;0.0867;0.0876;0.0886;0.0896;0.0906;0.0917;0.0927;0.0937;0.0948;0.0958;0.0969;0.0979;0.099;0.1;0.1011;0.1021;0.1032;0.1042;0.1053;0.1064;0.1074;0.1085;0.1095;0.1106;0.1117;0.1128;0.1139;0.115;0.116;0.1171;0.1182;0.1193;0.1204;0.1215;0.1225;0.1236;0.1246;0.1257;0.1267;0.1277;0.1287;0.1297;0.1308;0.1318;0.1328;0.1338;0.1348;0.1358;0.1368;0.1379;0.1389;0.1399;0.141;0.142;0.143;0.144;0.145;0.146;0.147;0.148;0.1489;0.1499;0.1508;0.1517;0.1525;0.1534;0.1543;0.1551;0.1559;0.1567;0.1575;0.1583;0.1591;0.1599;0.1606;0.1614;0.1621;0.1628;0.1635;0.1642;0.1649;0.1656;0.1663;0.167;0.1676;0.1683;0.1689;0.1696;0.1703;0.1709;0.1716;0.1723;0.1729;0.1736;0.1743;0.175;0.1757;0.1763;0.177;0.1777;0.1783;0.1789;0.1796;0.1802;0.1808;0.1814;0.182;0.1826;0.1831;0.1837;0.1843;0.1848;0.1854;0.1859;0.1865;0.1871;0.1876;0.1882;0.1887;0.1893;0.1899;0.1904;0.191;0.1915;0.1921;0.1926;0.1932;0.1937;0.1942;0.1948;0.1953;0.1958;0.1963;0.1969;0.1974;0.1979;0.1985;0.199;0.1995;0.2001;0.2006;0.2012;0.2017;0.2023;0.2028;0.2033;0.2038;0.2043;0.2048;0.2053;0.2057;0.2062;0.2066;0.207;0.2074;0.2078;0.2082;0.2086;0.209;0.2093;0.2097;0.21;0.2103;0.2105;0.2108;0.211;0.2112;0.2114;0.2115;0.2116;0.2117;0.2117;0.2118;0.2118;0.2117;0.2117;0.2116;0.2116;0.2115;0.2114;0.2113;0.2111;0.2109;0.2108;0.2105;0.2103;0.21;0.2096;0.2093;0.2089;0.2084;0.2079;0.2074;0.2069;0.2063;0.2057;0.2051;0.2045;0.2039;0.2032;0.2025;0.2019;0.2012;0.2005;0.1998;0.1992;0.1985;0.1978;0.1971;0.1964;0.1957;0.1951;0.1944;0.1937;0.1931;0.1924;0.1918;0.1912;0.1905;0.1899;0.1894;0.1888;0.1882;0.1877;0.1871;0.1866;0.1861;0.1856;0.1851;0.1846;0.1842;0.1837;0.1833;0.1829;0.1825;0.1821;0.1818;0.1815;0.1812;0.181;0.1808;0.1806;0.1805;0.1804;0.1803;0.1802;0.1801;0.18;0.1799;0.1797;0.1796;0.1794;0.1792;0.179;0.1787;0.1784;0.1781;0.1777;0.1773;0.1769;0.1765;0.176;0.1756;0.1751;0.1747;0.1742;0.1737;0.1733;0.1729;0.1724;0.172;0.1716;0.1712;0.1708;0.1704;0.17;0.1697;0.1694;0.169;0.1687;0.1684;0.168;0.1674;0.1666;0.1656;0.1644;0.1631;0.1621;0.1616;0.1614;0.1615;0.1617;0.1621;0.1625;0.1629;0.1632;0.1636;0.1639;0.1642;0.1645;0.1647;0.165;0.1653;0.1655;0.1658;0.166;0.1662;0.1665;0.1667;0.1669;0.1671;0.1673;0.1675;0.1677;0.1679;0.1681;0.1683;0.1685;0.1687;0.1689;0.1691;0.1693;0.1695;0.1697;0.1699;0.1701;0.1702;0.1704;0.1706;0.1707;0.1709;0.171;0.1711;0.1712;0.1713;0.1714;0.1714;0.1714;0.1714;0.1713;0.1712;0.171;0.1708;0.1706;0.1703;0.1699;0.1695;0.169;0.1684;0.1678;0.1671;0.1663;0.1655;0.1645;0.1635;0.1624;0.1612;0.1599;0.1584;0.1569;0.1553;0.1536;0.1517;0.1498;0.1477;0.1455;0.1433;0.1409;0.1384;0.1358;0.1331;0.1303;0.1274;0.1244;0.1213;0.1181;0.1149;0.1116;0.1082;0.1047;0.1013;0.0977;0.0942;0.0906;0.087;0.0834;0.0797;0.0761;0.0726;0.069;0.0655;0.0621;0.0587;0.0553;0.0521;0.0489;0.0458;0.0428;0.0399;0.0371;0.0344;0.0318;0.0294;0.0271;0.0248;0.0228;0.0208;0.019;0.0173;0.0157;0.0143;0.0129;0.0117;0.0106;0.0096;0.0087;0.0079;0.0072;0.0066;0.006;0.0055;0.0051;0.0048;0.0045;0.0043;0.0041;0.004;0.0039;0.0038;0.0038;0.0038;0.0038;0.0039;0.004;0.0041;0.0042;0.0043;0.0045;0.0046;0.0048;0.0049;0.0051;0.0053;0.0054;0.0056;0.0058;0.0059;0.0061;0.0063;0.0064;0.0066;0.0067;0.0068;0.007;0.0071;0.0072;0.0073;0.0074;0.0075;0.0075;0.0076;0.0077;0.0077;0.0078;0.0078;0.0079;0.0079;0.0079;0.0079;0.008;0.008;0.008;0.008;0.008;0.0081;0.0081;0.0081;0.0081;0.0081;0.0081;0.0081;0.0082;0.0082;0.0082;0.0083;0.0083;0.0084;0.0085;0.0086;0.0087;0.0089;0.009;0.0092;0.0094;0.0096;0.0099;0.0101;0.0104;0.0107;0.011;0.0113;0.0116;0.0119;0.0123;0.0126;0.013;0.0134;0.0138;0.0142;0.0146;0.015;0.0154;0.0159;0.0163;0.0168;0.0173;0.0177;0.0182;0.0187;0.0192;0.0197;0.0202;0.0207;0.0212;0.0216;0.0221;0.0226;0.0231;0.0235;0.024;0.0245;0.0249;0.0254;0.0258;0.0262;0.0267;0.0271;0.0275;0.028;0.0284;0.0288;0.0292;0.0296;0.03;0.0304;0.0308;0.0313;0.0317;0.0321;0.0325;0.0329;0.0333;0.0337;0.0342;0.0346;0.035;0.0354;0.0358;0.0362;0.0367;0.0371;0.0375;0.0379;0.0383;0.0387;0.0392;0.0396;0.04;0.0404;0.0408;0.0413;0.0417;0.0421;0.0425;0.0429;0.0434;0.0438;0.0442;0.0446;0.045;0.0454;0.0458;0.0462;0.0466;0.047;0.0474;0.0477;0.0481;0.0485;0.0488;0.0492;0.0495;0.0498;0.0501;0.0504;0.0507;0.051;0.0513;0.0516;0.0518;0.0521;0.0523;0.0526;0.0528;0.0531;0.0533;0.0536;0.0538;0.054;0.0542;0.0545;0.0547;0.0549;0.0551;0.0553;0.0555;0.0557;0.0559;0.0561;0.0563;0.0565;0.0567;0.057;0.0572;0.0575;0.0577;0.058;0.0583;0.0586;0.0589;0.0592;0.0596;0.06;0.0603;0.0607;0.0612;0.0616;0.062;0.0625;0.0629;0.0634;0.0638;0.0643;0.0648;0.0653;0.0658;0.0663;0.0668;0.0673;0.0678;0.0683;0.0688;0.0693;0.0699;0.0704;0.0709;0.0715;0.072;0.0726;0.0731;0.0736;0.0741;0.0747;0.0752;0.0757;0.0762;0.0766;0.0771;0.0776;0.078;0.0785;0.0789;0.0793;0.0796;0.08;0.0804;0.0807;0.081;0.0813;0.0816;0.0819;0.0822;0.0824;0.0827;0.083;0.0832;0.0834;0.0837;0.0839;0.0841;0.0843;0.0845;0.0847;0.085;0.0852;0.0854;0.0856;0.0857;0.0859;0.0861;0.0863;0.0865;0.0866;0.0868;0.087;0.0871;0.0873;0.0874;0.0875;0.0876;0.0877;0.0878;0.0878;0.0879;0.0879;0.088;0.088;0.088;0.088;0.088;0.088;0.088;0.088;0.0879;0.0879;0.0878;0.0877;0.0876;0.0875;0.0874;0.0872;0.087;0.0868;0.0867;0.0864;0.0862;0.086;0.0858;0.0856;0.0854;0.0851;0.0849;0.0847;0.0845;0.0843;0.0841;0.0838;0.0836;0.0834;0.0831;0.0829;0.0826;0.0823;0.082;0.0817;0.0813;0.081;0.0806;0.0802;0.0798;0.0794;0.079;0.0787;0.0783;0.0779;0.0775;0.0772;0.0768;0.0765;0.0761;0.0758;0.0755;0.0751;0.0748;0.0745;0.0742;0.0739;0.0735;0.0732;0.0728;0.0725;0.0721;0.0717;0.0713;0.0709;0.0704;0.07;0.0695;0.0691;0.0686;0.068;0.0675;0.067;0.0664;0.0658;0.0652;0.0645;0.0639;0.0632;0.0626;0.0619;0.0611;0.0604;0.0597;0.059;0.0582;0.0575;0.0568;0.0561;0.0554;0.0547;0.0541;0.0534;0.0528;0.0522;0.0516;0.0511;0.0505;0.05;0.0495;0.049;0.0486;0.0481;0.0477;0.0473;0.0469;0.0465;0.0462;0.0458;0.0454;0.0451;0.0448;0.0444;0.0441;0.0438;0.0434;0.043;0.0427;0.0423;0.0419;0.0415;0.0411;0.0407;0.0403;0.0399;0.0395;0.0391;0.0387;0.0383;0.0379;0.0376;0.0373;0.037;0.0368;0.0365;0.0363;0.0362;0.036;0.0359;0.0358;0.0357;0.0357;0.0356;0.0356;0.0355;0.0355;0.0355;0.0354;0.0354;0.0354;0.0353;0.0353;0.0353;0.0353;0.0353;0.0353;0.0353;0.0353;0.0353;0.0353;0.0354;0.0354;0.0354;0.0355;0.0355;0.0355;0.0356;0.0356;0.0357;0.0357;0.0358;0.0359;0.036;0.0361;0.0363;0.0365;0.0367;0.037;0.0373;0.0376;0.038;0.0384;0.0388;0.0392;0.0396;0.0401;0.0405;0.041;0.0415;0.0419;0.0424];           % Tabletran, 2001 coefficients (400 nm to 2400 nm) 


refl=Tablerefl(lambda-399); 
tran=Tabletran(lambda-399);   

