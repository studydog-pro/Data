si = [];
si.lamuda_CC=0.269;%碳捕集能耗Mwh/t
si.eta_cc=0.85;%碳捕集效率
si.eta_H2=0.75;
si.eta=39.65;%kwh/kg
si.eta_MR=0.75;%MR效率
si.aerfa=0.1;%掺氢比10%
si.daita=2;%甲烷和氢气转换系数
si.rou_H2=39.062;%氢气密度（70MPa，39.062kg/m3）
si.rou_CH4=312.496;%甲烷密度
si.rou_mix=285.1526;%混合密度（10%）
si.M_mix=14.436;%以掺氢比（10%）的混合摩尔质量
si.lamuda_CH4=0.8;%甲烷燃烧反应生成物计量数
si.lamuda_H2=0.9;%氢气燃烧反应生成物计量数 
si.fai_CH4=0.888;%甲烷摩尔分数转换系数
si.fai_H2=0.112;%氢气摩尔分数转换系数
si.T_in=640;%HT进气口燃烧温度K
si.T_comCH4=2223;%甲烷稳态燃烧温度K
si.T_comH2=2293;%氢气稳态燃烧温度K
si.HeatCH4=39.8;%甲烷热值
si.HeatH2=12.75;%氢气热值
si.Linepack0=100;%管存初值
si.aerfa_nv_CH4=1;%甲烷体积换算系数
si.aerfa_nv_H2=1;%氢气体积换算系数
si.N_in_CH4=0.9;%甲烷初始摩尔分数
si.N_in_H2=0.1;%氢气初始摩尔分数
si.a=0.9;%燃煤机组单位出力产生二氧化碳量1t/mw
si.daita_e=0.2;%碳配额系数t/MW
si.lamuda=0.2;%碳补偿系数
si.beta=39.3;%碳交易价格
si.l=20;%碳交易区间长度20t
si.e=0.25;%碳交易价格增长幅度
si.Horizon=24;

si.M_air=28.9634;si.M_H2=2.01588;si.M_CH4=16.043;%摩尔质量
si.absorb=0.95;
si.tao=0.03;%GC 耗气率
si.fp=3.5;%天然气燃料成本系数3.5     CCS:1.7  MR:0.3
si.HTcost=1.5;%CHP的燃料成本系数  1.5
si.MRcost=2.9;%MR的成本系数1:0.15    2.783  2.789   2.66/2.66/2.9
si.CCScost=46;%CCS燃料成本系数154元/MW
si.NumEn=30;si.NumEl=41;
si.NumGn=20;si.NumGl=19;si.NumGC=3;si.NumGW=4;
si.NumHT=2;%HGMT数
si.NumPG=3;
si.inonoff=[1 1 0]';%initial state of thermal units
si.minup=[4 4 4];
si.mindown=si.minup;
si.c=[0 0 0];si.stcost=[2800 2800 2800];si.sdcost=[3200 3200 3200];
si.HTlocat_E=[22;25];si.HTlocat_G=[2;7];%HT在电网和气网的安装位置
si.PGlocat=[1;2;13];%nonNGU在电网的安装位置
si.data=xlsread('Power.xlsx','Line data','A2:D42');
PPD=xlsread('Power','Units data','A2:F4');%nonNGU的数据
LNO=xlsread('Power','Bus load','A2:B21');
GLD=xlsread('Gas.xlsx','Gas load','A2:D21');
GWD=xlsread('Gas.xlsx','Gas well','A2:E5');
si.TPD=xlsread('Gas.xlsx','Pipeline','A2:E17');
GCD=xlsread('Gas.xlsx','Compressor','A2:F4');

%氢混燃机 = 数量   电网位置  气网位置  最低耗气量  最高耗气量

HT=[
1	22	2	500	1500
2	25	7	500	1500
];
si.NumHT=2;
si.HTlocat_E=HT(:,2);si.HTlocat_G=HT(:,3);
si.VHTmin=zeros(si.NumGn,1);si.VHTmin(si.HTlocat_G)=HT(:,4);
si.VHTmax=zeros(si.NumGn,1);si.VHTmax(si.HTlocat_G)=HT(:,5);
si.PHTmin=zeros(si.NumEn,1);si.PHTmin(si.HTlocat_E)=0;
si.PHTmax=zeros(si.NumEn,1);si.PHTmax(si.HTlocat_E)=600;%HT的数量，在电网和气网中的位置,最大耗气量,最大最小功率（取值足够大）

%Number 	Bus at the power system 	Node at the gas system 	Conversion coefficient (kcf/MW)	P2G min	P2G max
P2G=[
1	10	4	1.5000	0	100
2	23	10	1.2000	0	100
];
si.NumP2G=2;
si.P2Glocat_E=P2G(:,2);si.P2Glocat_G=P2G(:,3);
si.P2Gmin=zeros(si.NumEn,1);si.P2Gmin(si.P2Glocat_E)=P2G(:,5);
si.P2Gmax=zeros(si.NumEn,1);si.P2Gmax(si.P2Glocat_E)=P2G(:,6);
si.P2G_Gmax=zeros(si.NumGn,1);si.P2G_Gmax(si.P2Glocat_G)=P2G(:,6);%P2G的数量，在电网和气网中的位置,最大功率

%Number 	Bus at the power system 	Node at the gas system 	WMR min	WMR max
MR=[
1	10	4	0	6000   0    100
2	23	10	0	6000   0    100
];
si.NumMR=2;
si.MRlocat_E=MR(:,2);si.MRlocat_G=MR(:,3);
si.GHMRmin=zeros(si.NumGn,1);si.GHMRmin(si.MRlocat_G)=MR(:,4);
si.GHMRmax=zeros(si.NumGn,1);si.GHMRmax(si.MRlocat_G)=MR(:,5);
si.PMRmin=zeros(si.NumEn,1);si.PMRmin(si.MRlocat_E)=MR(:,6);
si.PMRmax=zeros(si.NumEn,1);si.PMRmax(si.MRlocat_E)=MR(:,7);%MR的数量，在电网和气网中的位置，最大耗气量

%Number 	Bus at the power system 	CCS min	  CCS max
CCS=[
1	1	0	100
2	2	0	100
3	13	0	100
4	22	0	100
5	25	0	100 
];
si.NumCCS=5;si.CCSlocat=CCS(:,2);
si.CCSmin=zeros(si.NumEn,1);
si.CCSmax=zeros(si.NumEn,1);si.CCSmax(si.CCSlocat)=CCS(:,4);%CCS的数量，在电网中的位置,最大功率

%Number 	Bus at the power system 	V min	  V max
rich=[
1	1	0	10000  
2	2	0	10000
3	13	0	10000
4	22	0	10000
5	25	0	10000
];
si.Numrich=5;si.richlocat=rich(:,2);
si.richmin=zeros(si.NumEn,1);si.richmin(si.richlocat)=rich(:,3);
si.richmax=zeros(si.NumEn,1);si.richmax(si.richlocat)=rich(:,4);

poor=[
1	1	0	10000
2	2	0	10000
3	13	0	10000
4	22  0	10000
5	25	0	10000
];
si.Numpoor=5;si.poorlocat=poor(:,2);
si.poormin=zeros(si.NumEn,1);si.poormin(si.richlocat)=poor(:,3);
si.poormax=zeros(si.NumEn,1);si.poormax(si.poorlocat)=poor(:,4);
si.co2=0.1;%富液罐co2浓度

% 储氢罐=数量     位置      最小容量     最大容量（m3）
H_storage=[
    1    4  4000  20000 0 5000; 
    2    10  4000  20000 0 5000; 
    ];
 si.NumHys=2;
 si.Hyslocat=H_storage(:,2);
 si.Hysmin=zeros(si.NumGn,1); si.Hysmin(si.Hyslocat)=H_storage(:,3);
 si.Hysmax=zeros(si.NumGn,1); si.Hysmax(si.Hyslocat)=H_storage(:,4);
 si.GHysmin=zeros(si.NumGn,1); si.GHysmin(si.Hyslocat)=H_storage(:,5);
 si.GHysmax=zeros(si.NumGn,1); si.GHysmax(si.Hyslocat)=H_storage(:,6);%Hys的数量，位置，最大，最小储气量   %%（最大最小流量）
%% 

si.PGmin=zeros(si.NumPG,1);si.PGmin=PPD(1:3,6);%nonNGU最大最小出力
si.PGmax=zeros(si.NumPG,1);si.PGmax=PPD(1:3,5);
si.NGc=PPD(1:3,2);si.NGb=PPD(1:3,3);si.NGa=PPD(1:3,4);%nonNGU发电成本系数
si.ud=si.PGmax*0.5;%nonNGU爬坡速率

si.GWlocat=GWD(:,2);%气源安装位置
si.Gflowmax=si.TPD(:,5);%气支路潮流限制
si.pimin=GLD(:,4);si.pimax=GLD(:,3);%气节点压力最大最小
si.GCmax=GCD(:,6);si.C_in=GCD(:,2);si.C_out=GCD(:,3);si.C_ratiomin=GCD(:,4);si.C_ratiomax=GCD(:,5);%压缩机相关参数
si.GWmin=zeros(si.NumGn,1);si.GWmin(GWD(:,2))=0;%气源最大最小出力
si.GWmax=zeros(si.NumGn,1);si.GWmax(GWD(:,2))=GWD(:,3);
%% 
% 风机最小出力  风机额定容量（没用到）  风机所在电网的节点位置
PWDATA=[
    0    300 10; 
    0    300 23; 
    ];%弃风弃光 数量太多
% 光伏最小出力  光伏额定容量（没用到）  光伏所在电网的节点位置
PVDATA=[
    0    300 10; 
    0    300 23; 
    ];
si.PWlocat=PWDATA(:,3);%风机所在节点位置
si.PVlocat=PVDATA(:,3);%光伏所在节点位置

PW=zeros(si.NumEn,1);
PWrate=[0.67 0.68 0.80 0.82 0.82 0.92 0.98 1.00 0.88 0.86 0.84 0.84 0.82 0.79 0.82 0.85 0.87 0.89 0.85 0.80 0.89 0.85 0.87 0.88];
% PWrate=[0.67 0.68 0.80 0.82 0.82 0.92 0.98 1.00 0.88 0.86 0.84 0.84 0.82 0.79 0.82 0.85 0.87 0.89 0.85 0.80 0.89 0.85 0.87 0.88];%风电典型出力曲线
% PWrate=[0.67 0.68 0.80 0.82 0.82 0.92 0.98 1.00 0.88 0.86 0.84 0.84 0.84 0.84 0.86 0.88 1.00 0.98 0.92 0.82 0.82 0.80 0.68 0.67];
for t=1:si.Horizon
    PW(si.PWlocat,t)=250*PWrate(t);%实际风电功率曲线
end
si.PW=PW;

PV=zeros(si.NumEn,1);
PVrate=[0 0 0 0 0 0.1 0.3 0.55 0.65 0.85 0.99 1 0.9 0.8 0.6 0.5 0.3 0 0 0 0 0 0 0];
% PVrate=[0 0 0 0 0 0.1 0.21 0.33 0.45 0.61 0.65 0.65 0.63 0.55 0.49 0.45 0.33 0.21 0.1 0 0 0 0 0];%光伏典型出力曲线
for t=1:si.Horizon
    PV(si.PVlocat,t)=250*PVrate(t);%实际光伏功率曲线
end
si.PV=PV;

%% P=B*theta 直流潮流方程
si.Y=zeros(si.NumEn);
for i=1:si.NumEl
    p=si.data(i,1);q=si.data(i,2);si.Y(p,q)=si.Y(p,q)-1./si.data(i,3);si.Y(q,p)=si.Y(p,q);si.Y(q,q)=si.Y(q,q)+1./si.data(i,3);si.Y(p,p)=si.Y(p,p)+1./si.data(i,3);
end
si.M=zeros(si.NumEl,si.NumEn);
for i=1:si.NumEl
    p=si.data(i,1);q=si.data(i,2);si.M(i,p)=1;si.M(i,q)=-1;
end

%% 电气负荷曲线
% PL=[0.68 0.64 0.62 0.60 0.61 0.63 0.68 0.70 0.73 0.81 0.89 0.92 0.95 0.95 0.97 0.95 0.93 0.91 0.89 0.90 0.93 0.91 0.86 0.86]*1600;%电负荷典型曲线 1900
PL=[0.68 0.64 0.62 0.60 0.61 0.63 0.68 0.70 0.73 0.81 0.89 0.92 0.92 0.89 0.81 0.73 0.70 0.68 0.63 0.61 0.60 0.62 0.64 0.68]*1600;%电负荷典型曲线 1900
si.PL=zeros(si.NumEn,si.Horizon);
for t=1:si.Horizon
    for i=1:size(LNO,1)
        si.PL(LNO(i,1),t)=PL(t)*LNO(i,2);%实际电负荷曲线
    end
end
GLrate=[1.0 1.0 0.9 0.9 0.85 1.0 1.0 0.6 0.6 0.6 0.8 0.7 1.0 1.0 0.7 0.7 1.0 1.0 1.0 0.7 0.7 0.7 0.8 1.0];
% GLrate=[1.0 1.0 0.9 0.9 1.0 1.0 1.0 0.7 0.7 0.7 0.8 0.7 1.0 1.0 0.9 0.9 1.0 1.0 1.0 0.7 0.7 0.7 0.8 1.0];%气负荷典型曲线
% GLrate=[1.0 1.0 0.9 0.9 1.0 1.0 1.0 0.7	0.7	0.7	0.8	0.7	0.7 0.8 0.7 0.7 0.7 1.0 1.0 1.0 0.9 0.9 1.0 1.0];
GL=zeros(si.NumGn,1);
GL(GLD(:,1))=GLD(:,2)*0.5;
for t=1:si.Horizon
    si.GL(:,t)=GL'*GLrate(t);%实际气负荷曲线
end