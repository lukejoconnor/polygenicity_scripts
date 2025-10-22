% Match the target size estimates from Simons et al. bioRxiv with the
% estimated polygenicity

matches = [1 51; 3 56; 5 82; 6 1; 9 57; 10 70; ...
    11 69; 13 38; 15 4; 17 25; 20 19; 33 45; 35 24];
our_traits = {'Body mass index'
'Years of education'
'FVC (corrected for smoking)'
'Neuroticism'
'Blood pressure - diastolic'
'Height'
'Chronotype (morning person)'
'Age at menarche'
'Total protein level in blood'
'Waist-hip ratio (corrected for BMI)'
'FEV1/FVC (corrected for smoking)'
'Reaction time'
'Creatinine level'
'Schizophrenia'
'Platelet count'
'Sleep duration'
'IGF1 level in blood'
'General risk tolerance'
'Aspartate aminotransferase level'
'Red blood cell count'
'Insomnia'
'Atrial fibrillation'
'Hypothyroidism'
'Allergy or eczema'
'Reported drinks per week'
'Phosphate levels'
'Breast cancer'
'Inflammatory bowel disease'
'Age at menopause'
'Bipolar disorder'
'Alzheimer''s disease'
'Coronary artery disease'
'Triglycerides'
'LDL'
'HbA1c'
'Hair color'};

other_traits = {'Standing height'
'Mean platelet (thrombocyte) volume'
'Sitting height'
'Platelet count'
'Mean corpuscular volume'
'Mean corpuscular haemoglobin'
'Platelet distribution width'
'Platelet crit'
'Trunk fat-free mass'
'Trunk predicted mass'
'Whole body fat-free mass'
'Whole body water mass'
'Alkaline phosphatase (quantile)'
'Monocyte count'
'Basal metabolic rate'
'Red blood cell (erythrocyte) distribution width'
'Eosinophill percentage'
'Mean reticulocyte volume'
'Red blood cell (erythrocyte) count'
'Mean sphered cell volume'
'Leg fat-free mass (right)'
'Leg predicted mass (right)'
'Monocyte percentage'
'Glycated haemoglobin (quantile)'
'IGF-1 (quantile)'
'Arm fat-free mass (right)'
'Arm predicted mass (right)'
'High light scatter reticulocyte percentage'
'High light scatter reticulocyte count'
'Lymphocyte count'
'Reticulocyte percentage'
'Reticulocyte count'
'Impedance of whole body'
'Weight'
'White blood cell (leukocyte) count'
'Gamma glutamyltransferase (quantile)'
'Weight'
'Creatinine (quantile)'
'Heel quantitative ultrasound index (QUI), direct entry'
'Heel bone mineral density (BMD) T-score, automated'
'Heel bone mineral density (BMD)'
'Impedance of leg (left)'
'Haemoglobin concentration'
'SHBG (quantile)'
'Triglycerides (quantile)'
'Impedance of arm (left)'
'Neutrophill count'
'Heel Broadband ultrasound attenuation, direct entry'
'Lymphocyte percentage'
'Haematocrit percentage'
'Body mass index (BMI)'
'Body mass index (BMI)'
'Immature reticulocyte fraction'
'Trunk fat mass'
'Whole body fat mass'
'Forced vital capacity (FVC)'
'Total protein (quantile)'
'Neutrophill percentage'
'Hip circumference'
'Arm fat mass (left)'
'Leg fat mass (left)'
'Ankle spacing width'
'Forced vital capacity (FVC), Best measure'
'Body fat percentage'
'Forced expiratory volume in 1-second (FEV1), predicted'
'Arm fat percentage (left)'
'Trunk fat percentage'
'Leg fat percentage (left)'
'Forced expiratory volume in 1-second (FEV1)'
'Waist circumference'
'Pulse rate, automated reading'
'Albumin (quantile)'
'Alanine aminotransferase (quantile)'
'Calcium (quantile)'
'Forced expiratory volume in 1-second (FEV1), Best measure'
'Heel bone mineral density (BMD) (right)'
'Heel quantitative ultrasound index (QUI), direct entry (right)'
'Heel bone mineral density (BMD) T-score, automated (right)'
'Systolic blood pressure, automated reading'
'Urea (quantile)'
'Heel broadband ultrasound attenuation (right)'
'Diastolic blood pressure, automated reading'
'Mean corpuscular haemoglobin concentration'
'Ankle spacing width (right)'
'3mm weak meridian (left)'
'3mm strong meridian (left)'
'Age first had sexual intercourse'
'Basophill percentage'
'6mm strong meridian (left)'
'Hand grip strength (right)'
'6mm weak meridian (left)'
'Corneal resistance factor (right)'
'Peak expiratory flow (PEF)'
'Spherical power (left)'
'Birth weight'};

target_size = log10([24187569.38
4541279.664
26382728.46
7136459.658
4276933.937
4623428.266
5078705.942
9891167.998
30913458.89
35234064.03
39056521.56
35633491.91
11922455.44
9740073.516
49411675.44
4120265.708
7037887.568
3436081.033
10922903.89
3977480.389
42154815.43
41830047.98
4554648.61
6027226.371
12734958.26
43836310.09
39203904.9
5408979.935
5488321.425
12032895.77
4142940.905
4100944.35
76564130.1
63628086.57
14809871.52
4868420.887
72333313.64
15779779.44
6574741.857
6574741.857
6303878.196
83437049.58
12397034.4
7961202.111
6654139.539
56667235.55
10919792.75
8780959.122
12977511.65
12416707.27
145064050.6
156658898.2
3818103.483
107642675.2
96404270.17
46506788.47
15222074.21
16293343.06
62245716.15
132648658.9
131341775.5
17001277.69
54499081.54
171072226.8
15084770.91
178630788.9
159458968.6
143852311.6
51794256.32
162959604.3
10301331.75
9374834.949
8838107.757
4352552.57
42092786.56
6422521.019
6227901.319
6227901.319
79416828.64
2559236.031
7189069.093
39194625.17
704445.49
17525949.14
8780269.739
14811052.35
929881181.2
2124511.436
9547813.167
392773584.3
9408429.033
4246068.644
38897084.06
5436455.076
20433722.13]);

target_size_bounds = log10([16899990.32	31877790.31
3572716.138	5757992.641
18301594.64	34696131.35
5459335.847	8739019.68
3213479.703	5644184.302
3396001.118	6475171.45
3851475.025	6759678.03
7014499.171	12623602.05
21803038.65	41536063.51
24487145.57	45899354.67
27135377.01	51320754.46
26000118.39	47935678.26
6830347.75	18461138.2
6918684.253	13534131.6
32063540.37	66301795.14
3047354.585	5313574.402
5400996.144	9779144.603
2302721.994	4534188.904
7646920.488	15309185.06
2917634.945	5216858.991
27665194.98	56562438.19
26542317.82	57523977.66
3206032.414	6860069.74
4175371.688	8435669.614
8172795.362	18049022.16
28983887.91	58284355.55
26763633.04	56285015.91
3982940.525	8254831.487
3946112.677	8443315.449
8204911.788	18050700.32
2811946.128	6173302.722
2745054.243	5816735.888
49146126.53	126221197.7
37864195.51	94595160.82
10079537.93	21851358.36
3194066.811	7206230.743
42598628.87	108565252.8
9151403.155	24214835.12
4542745.817	9277693.693
4542745.817	9277693.693
4486374.404	9014103.081
53437205.7	153867660.2
8441971.09	17834020.18
4729944.901	14171659.04
3922946.694	11724447.83
40860455.5	88785180.37
7134462.346	16396866.8
5529451.956	12725020.43
8559327.127	18721166.96
9099723.007	20007199.7
65190197	355812280.3
72041322.74	372368137.3
2787371.494	5916871.715
59914585.66	195990789.2
43562733.34	168470713.6
26037206.56	71285101.22
8360927.099	29204933.71
10329694.11	22655797.61
31748164.89	97527756.19
50594469.93	306328282.5
63830137.44	281772152
11270978.2	24093440.14
34011813.46	90404559.53
78691833.12	354553903.7
9703096.957	22200714.48
65267000.04	473733384.7
73105639.81	341121137.5
64422605.07	318958534.5
24248146.53	85894783.26
66635412.1	460449296.3
5769899.25	15922373.43
5532172.711	18516216.96
5476980.363	14337260.41
2612059.144	6610942.337
19883269.39	78518314.85
4114680.207	10721131.45
3898018.569	9611559.461
3898018.569	9611559.461
44613037.66	155965882.1
1708497.346	3751301.047
3489872.384	13342776.14
19960079.81	78112958.04
493757.643	1375823.428
7496715.304	31116734.84
4007179.015	21281327.83
5899994.044	34796050.28
439410578.7	2156125785
1257918.881	4499170.818
4431611.182	22505500.35
121701611.5	1061822597
4082188.148	20722961.26
1948132.801	8710726.649
12849620.39	107654612.3
3006587.2	10813679.69
7348074.489	80207876.85]);

polygenicity = [4.717155161
4.766821481
4.768843198
4.676895407
4.423654919
4.316627264
4.732653364
4.572527838
4.111995326
4.327440408
4.243572249
4.609981305
4.500569552
4.741280946
3.796843435
4.664294384
4.255120001
4.911216849
4.012554111
3.82172359
4.81935411
4.399875856
4.000155686
4.225144559
4.71383577
4.12591139
3.884193489
3.86815755
4.461370709
4.79040152
4.14208895
3.996669707
3.6110006
2.79977335
3.888654863
2.612633517];

polygenicity_se = [0.042864703
0.055153035
0.026629532
0.053029107
0.054359699
0.045267614
0.068926558
0.061399273
0.073895515
0.053708412
0.067834154
0.093347762
0.078163683
0.033424223
0.089398359
0.061765713
0.071771252
0.127296992
0.0682769
0.068746081
0.384648537
0.071441373
0.111346706
0.071114144
0.063182894
0.120064456
0.120125001
0.108275584
0.124369051
0.149945015
0.306704944
0.09724051
0.190928225
0.184390443
0.083807796
0.261064209];

eff_polygenicity = [3.584001574
3.944340773
3.835271184
3.970830127
3.70848576
3.793203036
3.979374944
3.344030415
3.590240307
3.429548696
3.493239683
4.136209837
3.387712507
4.116247823
3.161177117
4.028325828
3.533465034
4.117506081
3.556900937
3.246190376
3.966770985
3.255033248
3.214416808
3.582396599
4.259151121
2.959435776
3.089193357
2.721800566
3.553451021
3.953113911
2.902412646
3.040116692
2.461228339
1.950466311
2.862455866
1.845230432];

xneg = target_size_bounds(:,1)-target_size;
xpos = target_size_bounds(:,2)-target_size;



matching_traits(:,1) = our_traits(matches(:,1));
matching_traits(:,2) = other_traits(matches(:,2));

errorbar(target_size(matches(:,2)), polygenicity(matches(:,1)), ...
    2*polygenicity_se(matches(:,1)), 2*polygenicity_se(matches(:,1)), ...
    xneg(matches(:,2)), xpos(matches(:,2)), '.', 'CapSize', 0);
xlabel('Target size'); ylabel('Entropy polygenicity')
xlim([6.5 9]); ylim([3 5.5])
set(gca,'XTick',7:9,'XTickLabel',{'10^7','10^8','10^9'})
set(gca,'YTick',3:5,'YTickLabel',{'1000','10^4','10^5'})



