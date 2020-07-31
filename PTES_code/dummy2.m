clc
clear all

N=[0.372254335	0.162472338
0.372272113	0.162382629
0.372272131	0.161160054
0.372273407	0.161920612
0.37227406	0.163513905
0.372553051	0.161534066
0.372760386	0.160417784
0.373177809	0.161561564
0.373576936	0.161504022
0.373906423	0.161021705
0.374008065	0.160360026
0.374863861	0.160437037
0.375885078	0.158997267
0.376256797	0.160870032
0.37641553	0.160079591
0.376572791	0.159443486
0.376996177	0.158996368
0.378111928	0.158708241
0.378156194	0.157670094
0.37821014	0.16006064
0.378380859	0.159342558
0.379776136	0.159277365
0.379877986	0.158690462
0.382754497	0.15910153
0.382766627	0.158540564
0.384266249	0.158137514
0.384354322	0.157872772
0.38563946	0.157206046
0.386016348	0.157089008
0.386587664	0.155661817
0.386640827	0.157835234
0.386692245	0.156821425
0.386980646	0.157080279
0.387766533	0.15650284
0.39057031	0.155787549
0.391866442	0.156207289
0.393754938	0.155200004
0.395448309	0.155724907
0.395515039	0.155040155
0.395608904	0.155150867
0.395861314	0.153886784
0.395938025	0.154160747
0.396380682	0.155028603
0.397631862	0.152900072
0.397696759	0.153888211
0.397727275	0.153427613
0.398001945	0.151971649
0.398018267	0.154446273
0.398039138	0.153328127
0.40305853	0.153375254
]

P=[0.404838075	0.154431012
0.419309221	0.152721771
0.398213366	0.155920599
0.373069549	0.160969768
0.387361288	0.158278528
0.412544322	0.152969848
0.430392993	0.151244226
0.385237978	0.160335854
0.41145496	0.154041304
0.394335965	0.156577455
0.399047399	0.155489315
0.386805144	0.159665655
0.411673949	0.1536214
0.389898318	0.156844969
0.396704291	0.156391887
]

figure(1)
plot(N(:,1),N(:,2),'d');
xlim([0.36 0.44])
ylim([0.15 0.164])
grid on; xlabel('1- Roundtrip efficiency'); ylabel('LCOS [USD/ kWh]');
legend('NSGA II');

figure(2)
plot(P(:,1),P(:,2),'ok');
xlim([0.36 0.44])
ylim([0.15 0.164])
grid on; xlabel('1- Roundtrip efficiency'); ylabel('LCOS [USD/ kWh]');
legend('MOPSO');




