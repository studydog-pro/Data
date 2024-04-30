function  [mas]=solve_master(si,penaltyfactor,record_Gflow,record_fea_flag);
%% 变量定义
iter = length(record_fea_flag);
slack=sdpvar(si.NumGl,si.Horizon,'full');%松弛变量

U=sdpvar(si.NumEn,si.Horizon,'full');%节点电压变量
PDN=sdpvar(si.NumEn,si.NumEn,si.Horizon,'full');QDN=sdpvar(si.NumEn,si.NumEn,si.Horizon,'full');%支路有功无功变量
Pb=sdpvar(si.Horizon,1);Ps=sdpvar(si.Horizon,1);%从上级电网买电、卖电变量
ub=binvar(1,si.Horizon);us=binvar(1,si.Horizon);%上级购电、售电二进制变量
dPW=sdpvar(si.NumEn,si.Horizon,'full');%弃风变量

SVC=sdpvar(si.NumEn,si.Horizon,'full');
PW_g=sdpvar(si.NumEn,si.Horizon,'full');PV_g=sdpvar(si.NumEn,si.Horizon,'full');%风光上网功率 
PW_cp=sdpvar(si.NumEn,si.Horizon,'full');PV_cp=sdpvar(si.NumEn,si.Horizon,'full');%风光供给P2G的功率 
% PW=zeros(si.NumEn,si.Horizon);PW(si.PWlocat,:)=si.PW;%每个节点上的风机出力???;%风光出力
PW=sdpvar(si.NumEn,si.Horizon,'full');
PW_1=zeros(si.NumEn,si.Horizon);PW_1(si.PWlocat,:)=si.PW;
P_P2G=sdpvar(si.NumEn,si.Horizon,'full');%P2G耗电量 
P_P2G_G=sdpvar(si.NumGn,si.Horizon,'full'); %
PMR=sdpvar(si.NumEn,si.Horizon,'full');%MR设备耗电量
P_de=sdpvar(si.NumEn,si.Horizon,'full');%海水淡化装置
PG=sdpvar(si.NumEn,si.Horizon,'full');

pia=sdpvar(si.NumGn,si.Horizon,'full');%气节点气压 
GW=sdpvar(si.NumGn,si.Horizon,'full');%气源出气量
GC=sdpvar(si.NumGn,si.Horizon,'full');%压缩机 
Gflow=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路气流 
GHT=sdpvar(si.NumGn,si.Horizon,'full');%HT的耗气量
GP2G=sdpvar(si.NumGn,si.Horizon,'full');%P2G生产氢气
GHMR=sdpvar(si.NumGn,si.Horizon,'full');%MR消耗的氢气流量 
GMR=sdpvar(si.NumGn,si.Horizon,'full');%MR产生甲烷流量
GHNG=sdpvar(si.NumGn,si.Horizon,'full');%P2G向气网输送的氢气
GHys_in=sdpvar(si.NumGn,si.Horizon,'full');GHys_out=sdpvar(si.NumGn,si.Horizon,'full');%储氢罐的流入、流出流量 
SHys=sdpvar(si.NumGn,si.Horizon,'full');%储氢罐容量 
uHysch=binvar(si.NumGn,si.Horizon,'full');uHysdch=binvar(si.NumGn,si.Horizon,'full');%表征储氢罐充放的01状态标识
GL=si.GL;%气负荷

PHT_E=sdpvar(si.NumEn,si.Horizon,'full');%HT的有功（电网维度） 
QHT_E=sdpvar(si.NumEn,si.Horizon,'full');%HT的无功

Direction=sdpvar(1,1);%
gamma=sdpvar(1,1);%气网方向指示变量
Prs_square=sdpvar(si.NumGn,si.Horizon,'full');%压力平方
PHI=sdpvar(si.NumGl,si.Horizon,'full');%辅助变量 

Gflow_1=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Gflow_2=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路正向、反向气流
Gflow_in1=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Gflow_out1=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路流入、流出气流
Gflow_in2=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Gflow_out2=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路反向流入、流出气流
Linepack=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Linepack_D=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%管存

Ratio=sdpvar(1,si.Horizon);H_node=sdpvar(1,si.Horizon);%掺氢比相关变量
W1=sdpvar(si.NumGn,si.Horizon,'full');W2=sdpvar(si.NumGl,si.Horizon,'full');%GHT、Gflow*H_node辅助变量
W3=sdpvar(si.NumGn,si.Horizon,'full');%GHT*Ratio辅助变量
Ratio_square=sdpvar(1,si.Horizon);
W4=sdpvar(si.NumGn,si.Horizon,'full');%掺氢比平方  掺氢比平方*GHT辅助变量

temNS=sdpvar(si.Node,si.Horizon,'full');%节点温度 供应节点
temNR=sdpvar(si.Node,si.Horizon,'full');%节点温度 返回节点
temPSout=sdpvar(si.Pipe,si.Horizon,'full');%出口侧管道温度 供应管道
temPRout=sdpvar(si.Pipe,si.Horizon,'full');%出口侧管道温度 返回管道
temPSin=sdpvar(si.Pipe,si.Horizon,'full');%注入测管道温度 供应管道
temPRin=sdpvar(si.Pipe,si.Horizon,'full');%注入测管道温度 返回管道

HRB=sdpvar(si.NumEn,si.Horizon,'full');%余热回收
PEB=sdpvar(si.NumEn,si.Horizon,'full');HEB=sdpvar(si.NumEn,si.Horizon,'full');%电锅炉

ES=sdpvar(si.NumEn,si.Horizon,'full');Pch=sdpvar(si.NumEn,si.Horizon,'full');Pdch=sdpvar(si.NumEn,si.Horizon,'full');%蓄电池剩余电量、充电功率、放电功率变量
uEch=binvar(si.NumEn,si.Horizon,'full');uEdch=binvar(si.NumEn,si.Horizon,'full');%蓄电池是否接入，1接入0不接入

C_HT=sdpvar(si.NumGn,si.Horizon,'full');
C_CCS1=sdpvar(si.NumEn,si.Horizon,'full');C_CCS2=sdpvar(si.NumEn,si.Horizon,'full');
v_rich=sdpvar(si.NumEn,si.Horizon,'full');v_poor=sdpvar(si.NumEn,si.Horizon,'full');
V_rich=sdpvar(si.NumEn,si.Horizon,'full');V_poor=sdpvar(si.NumEn,si.Horizon,'full');
C_rich=sdpvar(si.NumEn,si.Horizon,'full');
C_MR=sdpvar(si.NumEn,si.Horizon,'full');
PCCS=sdpvar(si.NumEn,si.Horizon,'full');%碳捕集装置耗电量 
C_cost=sdpvar(1,1);%碳交易成本

onoff=binvar(si.NumP2G,si.Horizon,'full');%on-off state of units
shut=binvar(si.NumP2G,si.Horizon,'full');%shut-down state of units
start=binvar(si.NumP2G,si.Horizon,'full');%start-up state of units

B1=binvar(si.NumP2G,si.Horizon);B2=binvar(si.NumP2G,si.Horizon);B3=binvar(si.NumP2G,si.Horizon);
omiga1=sdpvar(si.NumP2G,si.Horizon);omiga2=sdpvar(si.NumP2G,si.Horizon);omiga3=sdpvar(si.NumP2G,si.Horizon);omiga4=sdpvar(si.NumP2G,si.Horizon);
eta_H2=sdpvar(si.NumP2G,si.Horizon,'full');%产氢效率

%% Electricity grid
Constraints = [];
Constraints = [Constraints, U(1,:) == 12.66];%参考节点的电压始终为1pu ?
% Constraints = [Constraints, -si.QTiemax<=QDN(1,2,:)<=si.QTiemax];% 配网与上级电网的无功功率交互限制?
Constraints = [Constraints, PDN(1,1,:)==0;QDN(1,1,:)==0];
Constraints = [Constraints, PDN(1,1,:)==PDN(1,2,:)];% Noted 不用管
Constraints = [Constraints, QDN(1,1,:)==QDN(1,2,:)];% Noted 不用管

% for t=1:si.Horizon
% Constraints = [Constraints, PDN(1,2,t)==Pb(t)-Ps(t)];% 配网与上级电网交互的有功功率
% Constraints = [Constraints, 0<=Pb(t)<=ub(t)*si.PTiemax]; % 配网与上级电网买电的功率交互限制 有功
% Constraints = [Constraints, 0<=Ps(t)<=us(t)*si.PTiemax]; % 配网与上级电网卖电的功率交互限制 有功
% Constraints = [Constraints, ub(t)+us(t)<=1]; %同一时刻，不能同时买电又卖电
% end

for t=1:si.Horizon
%     Constraints = [Constraints,PHT_E(si.HTlocat_E,t)==1.25e-4*(((1-Ratio(1,t))^2*(0.716*(1-Ratio(1,t))+0.089*Ratio(1,t))*Ratio_help(1,t))*35061+11340*(1-Ratio(1,t))+(Ratio(1,t)^2*(0.716*(1-Ratio(1,t))+0.089*Ratio(1,t))*Ratio_help(1,t))*3848+1340*Ratio(1,t))*GHT(si.HTlocat_G,t)];%HT发电功率进行转换（其余行进行清零）0.1485679853
%     Constraints = [Constraints,PHT_E(si.HTlocat_E,t)==1.25e-4*(((1-Ratio(1,t))^2)*1210.53+11340*(1-Ratio(1,t))+(Ratio(1,t)^2)*132.86+1340*Ratio(1,t))*GHT(si.HTlocat_G,t)];%HT发电功率进行转换（其余行进行清零）0.1485679853
%     Constraints = [Constraints,PHT_E(si.HTlocat_E,t)==2e-5*(12550.53*GHT(si.HTlocat_G,t)-12421.06*GHT(si.HTlocat_G,t)*Ratio(1,t)-1077.67*GHT(si.HTlocat_G,t)*Ratio_square(1,t))];
    Constraints = [Constraints,PHT_E(si.HTlocat_E,t)==2e-5*(12550.53*GHT(si.HTlocat_G,t)-12421.06*W3(si.HTlocat_G,t)-1077.67*W4(si.HTlocat_G,t))];
%     Constraints = [Constraints,PHT_E(si.HTlocat_E,t)==0.1485679853*GHT(si.HTlocat_G,t)];
    Constraints = [Constraints,PMR(si.MRlocat_E,t)==0.02*GMR(si.MRlocat_G,t)];%MR耗电
    Constraints = [Constraints,C_HT(si.HTlocat_G,t)==0.2*si.lamuda_CH4*si.aerfa_nv_CH4*si.N_in_CH4*GHT(si.HTlocat_G,t)];%HT碳排量
end

for t=1:si.Horizon
    for i=1:si.NumGn
        Constraints = [Constraints,W3(i,t) >= si.GHTmin*Ratio(1,t)+GHT(i,t)*si.Ratiomin-si.GHTmin*si.Ratiomin];%McCormick 处理GHT*Ratio 
        Constraints = [Constraints,W3(i,t) >= si.GHTmax*Ratio(1,t)+GHT(i,t)*si.Ratiomax-si.GHTmax*si.Ratiomax];
        Constraints = [Constraints,W3(i,t) <= si.GHTmax*Ratio(1,t)+GHT(i,t)*si.Ratiomin-si.GHTmax*si.Ratiomin];
        Constraints = [Constraints,W3(i,t) <= si.GHTmin*Ratio(1,t)+GHT(i,t)*si.Ratiomax-si.GHTmin*si.Ratiomax];

        Constraints = [Constraints,W4(i,t) >= si.GHTmin*Ratio_square(1,t)+GHT(i,t)*(si.Ratiomin^2)-si.GHTmin*(si.Ratiomin^2)];%McCormick 处理GHT*Ratio^2
        Constraints = [Constraints,W4(i,t) >= si.GHTmax*Ratio_square(1,t)+GHT(i,t)*(si.Ratiomax^2)-si.GHTmax*(si.Ratiomax^2)];
        Constraints = [Constraints,W4(i,t) <= si.GHTmax*Ratio_square(1,t)+GHT(i,t)*(si.Ratiomin^2)-si.GHTmax*(si.Ratiomin^2)];
        Constraints = [Constraints,W4(i,t) <= si.GHTmin*Ratio_square(1,t)+GHT(i,t)*(si.Ratiomax^2)-si.GHTmin*(si.Ratiomax^2)];
    end
end

for t=1:si.Horizon
%     Constraints = [Constraints,C_CCS1(si.PGlocat,t)==si.absorb*C_G(si.PGlocat,t)];
    Constraints = [Constraints,C_CCS1(si.HTlocat_E,t)==si.absorb*C_HT(si.HTlocat_G,t)];%吸收塔吸收的co2
    Constraints = [Constraints,C_CCS2(si.CCSlocat,t)==C_CCS1(si.CCSlocat,t)+C_rich(si.CCSlocat,t)];%再生塔再生的co2
    Constraints = [Constraints,C_rich(si.CCSlocat,t)==si.co2*v_rich(si.richlocat,t)];%富液罐流出的co2
    Constraints = [Constraints,v_rich(si.richlocat,t)+v_poor(si.poorlocat,t)==0];%流入等于流出
    Constraints = [Constraints,si.richmin<=V_rich(:,t)<=si.richmax];
    Constraints = [Constraints,si.poormin<=V_poor(:,t)<=si.poormax];
end
for t=2:si.Horizon
    Constraints = [Constraints,V_rich(si.richlocat,t)==V_rich(si.richlocat,t-1)-v_rich(si.richlocat,t)];
    Constraints = [Constraints,V_poor(si.poorlocat,t)==V_poor(si.poorlocat,t-1)-v_poor(si.poorlocat,t)];
end
v_rich(setdiff(1:si.NumEn,[si.richlocat]),:)=0;v_poor(setdiff(1:si.NumEn,[si.poorlocat]),:)=0;
Constraints = [Constraints,V_rich(:,1)==si.richmax*0.1-v_rich(:,1)];
Constraints = [Constraints,V_poor(:,1)==si.poormax*0.1-v_poor(:,1)];
Constraints = [Constraints,V_rich(:,si.Horizon)==si.richmax*0.1];
Constraints = [Constraints,V_poor(:,si.Horizon)==si.poormax*0.1];%贫、富液罐的容量限制

for t=1:si.Horizon
    Constraints = [Constraints,PCCS(si.HTlocat_E,t)==si.lamuda_CC*si.eta_cc*C_CCS2(si.CCSlocat,t)];%CCS在HT机组处出力
%     Constraints = [Constraints,PCCS(si.PGlocat,t)==si.lamuda_CC*si.eta_cc*C_CCS2(si.PGlocat,t)];%CCS在火电机组处出力
end


PHT_E(setdiff(1:si.NumEn,[si.HTlocat_E]),:)=0;
QHT_E(setdiff(1:si.NumEn,[si.HTlocat_E]),:)=0;%HT功率无关行清零

for t=1:si.Horizon
    Constraints = [Constraints, si.QHTmin <= QHT_E(si.HTlocat_E,t) <= si.QHTmax]; %微燃机无功出力限制
    Constraints = [Constraints, si.PEBmin <= PEB(:,t) <= si.PEBmax]; %电锅炉有功出力限制
    Constraints = [Constraints, si.Vmin <= U(:,t) <= si.Vmax];%节点电压限制

    Constraints = [Constraints,si.P2Gmin<=P_P2G(:,t)<=si.P2Gmax];
    Constraints = [Constraints,P_P2G(:,t)==PW_cp(:,t)];%P2G装置由风光供电
    Constraints = [Constraints,si.PMRmin<=PMR(:,t)<=si.PMRmax];%MR设备
    Constraints = [Constraints,0<=PW_g(:,t)];%上网功率
    Constraints = [Constraints,0<=PW_cp(:,t)];%给P2G供气功率
    Constraints = [Constraints,PW(:,t)==PW_g(:,t)+PW_cp(:,t)];%风电出力限制
    Constraints = [Constraints,0<=PW(:,t)<=PW_1(:,t)];
    Constraints = [Constraints,si.CCSmin<=PCCS(:,t)<=si.CCSmax];%碳捕集装置耗能上下限
    Constraints = [Constraints,P_de(:,t)==P_P2G(:,t)*si.eta_de];%海水淡化电量
    Constraints = [Constraints,si.SVCmin<=SVC(:,t)<=si.SVCmax];%SVC
%     Constraints = [Constraints,si.PGmin<=PG(:,t)<=si.PGmax];
end

for t=1:si.Horizon
    if t==1
        Constraints=[Constraints,start(:,t)-shut(:,t)==onoff(:,t)-si.inonoff];
        Constraints=[Constraints,start(:,t)+shut(:,t)<=ones(si.NumP2G,1)];
    else
        Constraints=[Constraints,start(:,t)-shut(:,t)==onoff(:,t)-onoff(:,t-1)];
        Constraints=[Constraints,start(:,t)+shut(:,t)<=ones(si.NumP2G,1)];
    end
end
for t=2:si.Horizon
    for unit=1:si.NumP2G
        indicator=onoff(unit,t)-onoff(unit,t-1);
        range=t:min(si.Horizon,t+si.minup(unit)-1);
        Constraints=[Constraints,onoff(unit,range)>=indicator];
        indicator=onoff(unit,t-1)-onoff(unit,t);
        range=t:min(si.Horizon,t+si.mindown(unit)-1);
        Constraints=[Constraints,onoff(unit,range)<=1-indicator];
    end
end
for t=1:si.Horizon
    Constraints=[Constraints,onoff(:,t).*si.P2G_min<=P_P2G(si.P2Glocat_E,t)<=onoff(:,t).*si.P2G_max];%EL功率限制（含启停）
end

P_P2G(setdiff(1:si.NumEn,[si.P2Glocat_E]),:)=0;

for t=1:si.Horizon
for i=1:si.NumEl
    ii=si.Network(i,1);%from_bus
    jj=si.Network(i,2);%to_bus
    rr=si.Network(i,3);%r
    xx=si.Network(i,4);%x
    if ii~=1
        Constraints = [Constraints, U(ii,t)==U(jj,t)+(PDN(ii,jj,t)*rr+QDN(ii,jj,t)*xx)/(12.66*1000)];%节点电压计算公式
    end
end
end

for t=1:si.Horizon
    for jj = 1:si.NumEn
        ii = si.rho(jj);
        Constraints = [Constraints,PDN(ii,jj,t)==si.PL(jj,t)-PHT_E(jj,t)-PW_g(jj,t)+PCCS(jj,t)+PMR(jj,t)+P_de(jj,t)+PEB(jj,t)+Pch(jj,t)-Pdch(jj,t)+(sum(si.A(jj,:).*PDN(jj,:,t)))];% 每个节点的有功功率平衡
        Constraints = [Constraints,QDN(ii,jj,t)==si.QL(jj,t)+SVC(jj,t)-QHT_E(jj,t)-(PW_g(jj,t))*si.K+sum(si.A(jj,:).*QDN(jj,:,t))];% 每个节点的无功功率平衡
    end
end


for t=1:si.Horizon
    Constraints = [Constraints,0<=Pch(:,t)<=si.PESmax.*uEch(:,t)];%蓄电池充电功率限制
    Constraints = [Constraints,0<=Pdch(:,t)<=si.PESmax.*uEdch(:,t)];%蓄电池放电功率限制
    Constraints = [Constraints,si.Emin<=ES(:,t)<=si.Emax];%蓄电池最大容量限制
    Constraints = [Constraints,uEch(:,t)+uEdch(:,t)<=1];%蓄电池充放电状态
end

for t=2:si.Horizon
    Constraints = [Constraints,ES(:,t)==ES(:,t-1)+Pch(:,t)*si.etaES-Pdch(:,t)/si.etaES];%蓄电池储存电量约束
end
 
Constraints = [Constraints,ES(:,1)==si.Emax*0.2+Pch(:,1)*si.etaES-Pdch(:,1)/si.etaES];%第一个时刻
Constraints = [Constraints,ES(:,si.Horizon)==si.Emax*0.2];%一日内蓄电池最后一个时刻的储存电量要等于今天刚开始时候的储存电量

%% GAS GRID include direction and linepack

for t=1:si.Horizon
    Constraints = [Constraints,si.Hysmin<=SHys(:,t)<=si.Hysmax];%储氢罐容量上下限
    Constraints = [Constraints,0<=GHys_in(:,t)<=si.GHysmax.*uHysch(:,t)];%储氢罐冲氢流量上下限
    Constraints = [Constraints,0<=GHys_out(:,t)<=si.GHysmax.*uHysdch(:,t)];%储氢罐放氢流量上下限
    Constraints = [Constraints,uHysch(:,t)+uHysdch(:,t)<=1];%储氢罐充放氢状态
end
for t=2:si.Horizon
    Constraints = [Constraints,SHys(:,t)==SHys(:,t-1)+GHys_in(:,t)*0.95-GHys_out(:,t)/0.95];%储气罐储存电量约束
end

Constraints = [Constraints,SHys(:,1)==si.Hysmax*0.2+GHys_in(:,1)*0.95-GHys_out(:,1)/0.95];%第一个时刻储氢罐容量
Constraints = [Constraints,SHys(:,si.Horizon)==si.Hysmax*0.2];%一日内储气罐最后一个时刻的储存热量要等于今天刚开始时候的储存气量

for t=1:si.Horizon
    Constraints = [Constraints,si.GHMRmin<=GHMR(:,t)<=si.GHMRmax];%MR进氢约束
    Constraints = [Constraints,si.VHTmin<=GHT(:,t)<=si.VHTmax];%氢燃混机燃烧气体体积约束
    Constraints = [Constraints,0<=GHNG(:,t)];%P2G注入气网的氢气流量必须大于0
end

for t=1:si.Horizon
    Constraints = [Constraints,si.pimin<=pia(:,t)<=si.pimax];%节点压力约束
end

for t=1:si.Horizon
    Constraints = [Constraints,H_node(1,t)==si.HeatCH4*(1-Ratio(1,t))+si.HeatH2*Ratio(1,t)];%节点热值
    Constraints = [Constraints,si.H_nodemin<=H_node(1,t)<=si.H_nodemax];
    Constraints = [Constraints,si.Ratiomin<=Ratio(1,t)<=si.Ratiomax];%掺氢比限制
    Constraints = [Constraints,Ratio_square(1,t)==Ratio(1,t)^2];%掺氢比平方
    Constraints = [Constraints,si.Ratiomin^2<=Ratio_square(1,t)<=si.Ratiomax^2];%掺氢比平方的限制
end

for t=1:si.Horizon
    Constraints = [Constraints,si.eta_H2min<=eta_H2(:,t)<=si.eta_H2max];%效率上下限
    Constraints = [Constraints,P_P2G(si.P2Glocat_E,t)./si.P2Gmax(si.P2Glocat_E)==0*omiga1(:,t)+(1/3)*omiga2(:,t)+(2/3)*omiga3(:,t)+1*omiga4(:,t)];%效率线性化
    Constraints = [Constraints,eta_H2(:,t)==0.8196*omiga1(:,t)+0.7527*omiga2(:,t)+0.6992*omiga3(:,t)+0.6590*omiga4(:,t)];
    Constraints = [Constraints,omiga1(:,t)+omiga2(:,t)+omiga3(:,t)+omiga4(:,t)==1];
    Constraints = [Constraints,B1(:,t)+B2(:,t)+B3(:,t)==1];
    Constraints = [Constraints,omiga1(:,t)<=B1(:,t)];
    Constraints = [Constraints,omiga2(:,t)<=B1(:,t)+B2(:,t)];
    Constraints = [Constraints,omiga3(:,t)<=B2(:,t)+B3(:,t)];
    Constraints = [Constraints,omiga4(:,t)<=B3(:,t)];
end


% d=binvar(2,si.Horizon);
dd=binvar(2,si.Horizon);
ddd=binvar(2,si.Horizon);%implise 中if elseif

for t=2:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        Constraints = [Constraints,Linepack(m,n,t)-Linepack(m,n,t-1)==Gflow_in1(m,n,t)-Gflow_out1(m,n,t)+Gflow_in2(m,n,t)-Gflow_out2(m,n,t)];%相邻断面管存
%         Constraints = [Constraints,Linepack_D(m,n,t)==Linepack(m,n,t)-Linepack(m,n,t-1)];
    end
end
for s=1:si.NumGl-si.NumGC
    m=si.TPD(s,2);
    n=si.TPD(s,3);
    Constraints = [Constraints,Linepack(m,n,1)-si.Linepack0==Gflow_in1(m,n,1)-Gflow_out1(m,n,1)+Gflow_in2(m,n,1)-Gflow_out2(m,n,1)];%相邻断面管存初值
%     Constraints = [Constraints,Linepack_D(m,n,1)==Linepack(m,n,1)-si.Linepack0];
end
for t=1:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        Constraints = [Constraints,Linepack(m,n,t)>=si.Linepack0];%管存约束
%         Constraints = [Constraints,Linepack_D(m,n,t)>=0];%是否必须大于0
    end
end

for t=1:si.Horizon
    Constraints = [Constraints,si.pimin.^2<=Prs_square(:,t)<=si.pimax.^2];
    Constraints = [Constraints,Prs_square(:,t)==pia(:,t).^2];%节点压力约束  平方
end

for t=1:si.Horizon
    Constraints = [Constraints,si.GWmin<=GW(:,t)<=si.GWmax];%气源出力限
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        W=si.TPD(s,4);
        Constraints = [Constraints,sum(dd(:,t))==1];
        Constraints = [Constraints,implies(dd(1,t),[pia(m,t) >= pia(n,t), gamma == 1])];
        Constraints = [Constraints,implies(dd(2,t),[pia(m,t) <= pia(n,t), gamma == -1])];%添加气体流动限制方向

        Constraints = [Constraints,0<=Gflow_1(m,n,t)<=gamma*si.Gflowmax(s)];%支路气流正向流量限制
        Constraints = [Constraints,0<=Gflow_2(m,n,t)<=(1-gamma)*si.Gflowmax(s)];%支路气流反向流量限制
        Constraints = [Constraints,0<=Gflow_in1(m,n,t)<=si.Gflowmax(s)];%支路气流流入流量限制m-n
        Constraints = [Constraints,0<=Gflow_out1(m,n,t)<=si.Gflowmax(s)];%支路气流流出流量限制m-n
        Constraints = [Constraints,0<=Gflow_in2(m,n,t)<=si.Gflowmax(s)];%支路气流流入流量限制n-m
        Constraints = [Constraints,0<=Gflow_out2(m,n,t)<=si.Gflowmax(s)];%支路气流流出流量限制n-m
        Constraints = [Constraints,Gflow_1(m,n,t)==(Gflow_in1(m ,n,t)+Gflow_out1(m,n,t))/2];
        Constraints = [Constraints,Gflow_2(m,n,t)==(Gflow_in2(m,n,t)+Gflow_out2(m,n,t))/2];
        Constraints = [Constraints,Gflow(m,n,t)==Gflow_1(m,n,t)-Gflow_2(m,n,t)];
        Constraints = [Constraints,-si.Gflowmax(s)<=Gflow(m,n,t)<=si.Gflowmax(s)];%支路气流流量限制

%         Constraints = [Constraints,norm([Gflow(m,n,t)/W;pia(n,t)],2)<=pia(m,t)];%Weymouth公式

        Constraints = [Constraints,PHI(s,t)>=Gflow(m,n,t)^2/W^2];% A
        Constraints = [Constraints,(1-gamma)/2*(si.pimin(m)^2-si.pimax(n)^2) <= Prs_square(m,t)-Prs_square(n,t)<=(1+gamma)/2*(si.pimax(m)^2-si.pimin(n)^2)];%E
        Constraints = [Constraints,PHI(s,t)>=Prs_square(n,t)-Prs_square(m,t)+(gamma+1)*(si.pimin(m)^2-si.pimax(n)^2)];
        Constraints = [Constraints,PHI(s,t)>=Prs_square(m,t)-Prs_square(n,t)+(gamma-1)*(si.pimax(m)^2-si.pimin(n)^2)];% B
        Constraints = [Constraints,PHI(s,t)<=Prs_square(n,t)-Prs_square(m,t)+(gamma+1)*(si.pimax(m)^2-si.pimin(n)^2)];
        Constraints = [Constraints,PHI(s,t)<=Prs_square(m,t)-Prs_square(n,t)+(gamma-1)*(si.pimin(m)^2-si.pimax(n)^2)];% C  节点压力 数值

        Constraints = [Constraints,Linepack(m,n,t)==0.6*(pia(m,t)+pia(n,t))/2];%管存 管存常数0.6

        Constraints = [Constraints,W2(s,t) >= 0*H_node(1,t)+Gflow(m,n,t)*si.H_nodemin-0*si.H_nodemin];%McCormick 处理Gflow*H_node
        Constraints = [Constraints,W2(s,t) >= si.Gflowmax(s)*H_node(1,t)+Gflow(m,n,t)*si.H_nodemax-si.Gflowmax(s)*si.H_nodemax];
        Constraints = [Constraints,W2(s,t) <= si.Gflowmax(s)*H_node(1,t)+Gflow(m,n,t)*si.H_nodemin-si.Gflowmax(s)*si.H_nodemin];
        Constraints = [Constraints,W2(s,t) <= 0*H_node(1,t)+Gflow(m,n,t)*si.H_nodemax-0*si.H_nodemax];

        if iter~=0
            Constraints = [Constraints,PHI(s,t)-(record_Gflow(m,n,t)^2+2*record_Gflow(m,n,t)*(Gflow(m,n,t)-record_Gflow(m,n,t)))/(W^2)<=slack(s,t)];%本质上 是限制 Gflow^2/W^2>=PHI 即slack变成0
            Constraints = [Constraints,slack(s,t)>=0];
        end
    end
    for i=1:si.NumGn
        if ismember(i,si.C_in)
            k=find(i==si.C_in);
            Constraints = [Constraints,si.C_ratiomin(k)*pia(si.C_in(k),t)<=pia(si.C_out(k),t)<=si.C_ratiomax(k)*pia(si.C_in(k),t)];%压缩机支路气压  A  (如果加上这个约束的话就相当于默认了gamma==1，既由输入节点流向输出节点)           
            Constraints = [Constraints,Gflow(si.C_in(k),si.C_out(k),t)==(1+si.tao)*GC(k,t)];%压缩机损耗
            Constraints = [Constraints,0<=GC(k,t)<=si.GCmax(k)];%压缩机出力限制
            Constraints = [Constraints,Linepack(si.C_in(k),si.C_out(k),t)==0];%压缩机管道管存
        end
    end
    GHT(setdiff(1:si.NumGn,[si.HTlocat_G]),:)=0;%非HT和EL节点的耗气量、产气量置为零
    GHNG(setdiff(1:si.NumGn,[si.P2Glocat_G]),:)=0;
    GMR(setdiff(1:si.NumGn,[si.MRlocat_G]),:)=0;
    for i=1:si.NumGn
        if ismember(i,si.P2Glocat_G)
            kk=find(i==si.P2Glocat_G);

            Constraints = [Constraints,GP2G(i,t)==(P_P2G(si.P2Glocat_E(kk),t)*eta_H2(kk,t)*3.6*500)/si.HeatH2];%EL电-气转换 10 4
            Constraints = [Constraints,GMR(i,t)==(4*si.eta_MR*GHMR(si.MRlocat_G(kk),t )*si.HeatH2)/si.HeatCH4];%MR产甲烷
            Constraints = [Constraints,GHNG(i,t)==GP2G(si.P2Glocat_G(kk),t)-GHMR(si.P2Glocat_G(kk),t)-GHys_in(si.P2Glocat_G(kk),t)];
        end
    end

    temp=sdpvar(si.NumGn,1);
    for i=1:si.NumGn
        Constraints = [Constraints,W1(i,t) >= si.GHTmin*H_node(1,t)+GHT(i,t)*si.H_nodemin-si.GHTmin*si.H_nodemin];%McCormick 处理GHT*H_node 
        Constraints = [Constraints,W1(i,t) >= si.GHTmax*H_node(1,t)+GHT(i,t)*si.H_nodemax-si.GHTmax*si.H_nodemax];
        Constraints = [Constraints,W1(i,t) <= si.GHTmax*H_node(1,t)+GHT(i,t)*si.H_nodemin-si.GHTmax*si.H_nodemin];
        Constraints = [Constraints,W1(i,t) <= si.GHTmin*H_node(1,t)+GHT(i,t)*si.H_nodemax-si.GHTmin*si.H_nodemax];
   
        temp(i)=GW(i,t)*si.HeatCH4+GMR(i,t)*si.HeatCH4+GHNG(i,t)*si.HeatH2+GHys_out(i,t)*si.HeatH2-si.GL(i,t)*H_node(1,t)-W1(i,t);%天然气节点净注入气流   维数
        for s=1:si.NumGl-si.NumGC
            m=si.TPD(s,2);
            n=si.TPD(s,3);
%             Constraints = [Constraints,sum(d(:,t))==1];
%             Constraints = [Constraints,implies(d(1,t),[pia(m,t) >= pia(n,t), Direction == 1])];
%             Constraints = [Constraints,implies(d(2,t),[pia(m,t) <= pia(n,t), Direction == 0])];%添加气体流动限制方向
            Direction = 1;
            if i==m
                temp(i)=temp(i)-Direction*W2(s,t)+(1-Direction)*W2(s,t);
%                 temp(i)=temp(i)-Gflow(m,n,t);
            elseif i==n
                temp(i)=temp(i)+Direction*W2(s,t)-(1-Direction)*W2(s,t);
%                 temp(i)=temp(i)+Gflow(m,n,t);
            end
        end
        if ismember(i,si.C_in)
            k=find(i==si.C_in);
            temp(i)=temp(i)-Gflow(si.C_in(k),si.C_out(k),t)*H_node(1,t);
        end
        if ismember(i,si.C_out)
            kk=find(i==si.C_out);
            temp(i)=temp(i)+Gflow(si.C_in(kk),si.C_out(kk),t)*H_node(1,t);
        end
        Constraints = [Constraints,temp(i)==0];%天然气节点气流平衡
    end
end

for t=1:si.Horizon
    Constraints = [Constraints,C_MR(si.MRlocat_E,t)==(GMR(si.MRlocat_G,t)/(1.397*1e-3))*44*5e-7];%MR消耗的CO2
    Constraints = [Constraints,0<=sum(C_MR(si.MRlocat_E,t))<=sum(si.eta_cc*C_CCS2(si.CCSlocat,t))];
    Constraints = [Constraints,sum(GHNG(si.P2Glocat_G,t))+sum(GHys_out(si.P2Glocat_G,t))==Ratio(1,t)*(sum(GHNG(si.P2Glocat_G,t))+sum(GHys_out(si.P2Glocat_G,t))+sum(GW(si.GWlocat,t))+sum(GMR(si.MRlocat_G,t)))];%控制掺氢比
end

C_STO=0;
for t=1:si.Horizon
    C_STO=C_STO+sum(si.eta_cc*C_CCS2(si.CCSlocat,t))-sum(C_MR(si.MRlocat_E,t));
end
 %% Heat grid
for t=1:si.Horizon
    Constraints = [Constraints,HRB(:,t)==PHT_E(:,t)*si.eta_re*si.eta_rb*si.beta_rb];
    Constraints = [Constraints,HEB(:,t)==PEB(:,t)*si.etaEB];
end
for t=1:si.Horizon
    Constraints=[Constraints,sum(HRB(si.HTlocat_E,t))+sum(HEB(si.EBlocat,t))==si.c*si.NFrate(1)*(temNS(1,t)-temNR(1,t))];
    Constraints=[Constraints,si.TSmin<=temNS(:,t)<=si.TSmax];
    Constraints=[Constraints,si.TRmin<=temNR(:,t)<=si.TRmax];
end
for i=8:si.Node %%only load node  从8开始
    for t=1:si.Horizon
        Constraints=[Constraints,si.HL(i,t)==si.c*si.NFrate(i)*(temNS(i,t)-temNR(i,t))];% Heating load node power balance
    end
end
for i=1:si.Pipe
    for t=1:si.Horizon
        Constraints=[Constraints,temPSout(i,t)-si.temp_H(t)==(temPSin(i,t)-si.temp_H(t))*(1-(si.lambda(i)*si.Length(i))/(si.c*si.SPFrate(i)))];
        Constraints=[Constraints,temPRout(i,t)-si.temp_H(t)==(temPRin(i,t)-si.temp_H(t))*(1-(si.lambda(i)*si.Length(i))/(si.c*si.RPFrate(i)))];
    end
end
for t=1:si.Horizon
    Constraints=[Constraints,temPRout(1,t)==temNR(1,t)];
    Constraints=[Constraints,temPRout(2,t)==temNR(2,t)];
    Constraints=[Constraints,temPRout(3,t)*si.RPFrate(3)+temPRout(7,t)*si.RPFrate(7)+temPRout(10,t)*si.RPFrate(10)==temNR(3,t)*(si.RPFrate(3)+si.RPFrate(7)+si.RPFrate(10))];
    Constraints=[Constraints,temPRout(4,t)*si.RPFrate(4)+temPRout(8,t)*si.RPFrate(8)+temPRout(11,t)*si.RPFrate(11)==temNR(4,t)*(si.RPFrate(4)+si.RPFrate(8)+si.RPFrate(11))];
    Constraints=[Constraints,temPRout(5,t)*si.RPFrate(5)+temPRout(9,t)*si.RPFrate(9)+temPRout(12,t)*si.RPFrate(12)==temNR(5,t)*(si.RPFrate(5)+si.RPFrate(9)+si.RPFrate(12))];
    Constraints=[Constraints,temPRout(6,t)*si.RPFrate(6)+temPRout(13,t)*si.RPFrate(13)==temNR(6,t)*(si.RPFrate(6)+si.RPFrate(13))];
    
    Constraints=[Constraints,temPSin(1,t)==temNS(1,t)];
    Constraints=[Constraints,temPSin(2,t)==temNS(2,t)];
    Constraints=[Constraints,temPSin(3,t)==temNS(3,t)];
    Constraints=[Constraints,temPSin(4,t)==temNS(4,t)];
    Constraints=[Constraints,temPSin(5,t)==temNS(5,t)];
    Constraints=[Constraints,temPSin(6,t)==temNS(6,t)];
    Constraints=[Constraints,temPSin(7,t)==temNS(3,t)];
    Constraints=[Constraints,temPSin(8,t)==temNS(4,t)];
    Constraints=[Constraints,temPSin(9,t)==temNS(5,t)];
    Constraints=[Constraints,temPSin(10,t)==temNS(3,t)];
    Constraints=[Constraints,temPSin(11,t)==temNS(4,t)];
    Constraints=[Constraints,temPSin(12,t)==temNS(5,t)];
    Constraints=[Constraints,temPSin(13,t)==temNS(6,t)];
end
for i=1:si.Pipe
    for t=1:si.Horizon
        Constraints=[Constraints,temPRin(i,t)==temNR(i+1,t)];
        Constraints=[Constraints,temPSout(i,t)==temNS(i+1,t)];
    end
end
%% Carbon trade
C_q=0;C_IES=0;C_a=0;C_CCS=0;
for t=1:si.Horizon
    C_q=C_q+si.daita_e*(sum(PHT_E(si.HTlocat_E,t)));%碳配额
end

for t=1:si.Horizon
    C_IES=C_IES+sum(C_HT(si.HTlocat_G,t))+sum((1-Ratio(1,t))*2.165*1e-3*si.GL(:,t));%全部碳排
end

for t=1:si.Horizon
    C_CCS=C_CCS+si.eta_cc*(sum(C_HT(si.HTlocat_G,t)));%碳捕集量
end
C_a=C_IES-C_CCS;%实际碳排

u=binvar(6,1);
Constraints = [Constraints,sum(u(:,1))==1];

Constraints = [Constraints,implies(u(1,1),[C_a <= C_q-si.l, C_cost==-si.beta*(1+2*si.lamuda)*(C_q-si.l-C_a)])];

Constraints = [Constraints,implies(u(2,1),[C_q-si.l <= C_a<=C_q, C_cost==-si.beta*(1+2*si.lamuda)*si.l-si.beta*(1+si.lamuda)*(C_q-C_a)])];

Constraints = [Constraints,implies(u(3,1),[C_q <= C_a <= C_q+si.l, C_cost==si.beta*(C_a-C_q)])];

Constraints = [Constraints,implies(u(4,1),[C_q+si.l <= C_a <= C_q+2*si.l, C_cost==si.beta*si.l+si.beta*(1+si.e)*(C_a-C_q-si.l)])];

Constraints = [Constraints,implies(u(5,1),[C_q+2*si.l <= C_a <= C_q+3*si.l, C_cost==si.beta*(2+si.e)*si.l+si.beta*(1+2*si.e)*(C_a-C_q-2*si.l)])];

Constraints = [Constraints,implies(u(6,1),[C_q+3*si.l <= C_a, C_cost==si.beta*(3+3*si.e)*si.l+si.beta*(1+3*si.e)*(C_a-C_q-3*si.l)])];


%% Objective function
obj1=0;obj2=0;obj3=0;obj4=0;obj5=0;obj6=0;obj7=0;obj8=0;obj9=0;obj9a=0;obj9b=0;
% for t=1:si.Horizon
% %     obj1=obj1+si.cb(t)*Pb(t)-si.cs(t)*Ps(t);%与上级电网购电和售电成本
% %       obj1=obj1+PG(si.PGlocat,t)'*diag(si.NGa)*PG(si.PGlocat,t)+si.NGb'*PG(si.PGlocat,t)+sum(si.NGc);
% end
for t=1:si.Horizon
    obj2=obj2+si.fp*sum(GW(si.GWlocat,t));%气源买气成本
end
for t=1:si.Horizon
    obj3=obj3+si.HTcost*sum(PHT_E(si.HTlocat_E,t));%HT运行成本
end

if iter==0
    for t=1:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        obj4=obj4+0.1*PHI(s,t);%weymouth等式松弛惩罚项   0.1取得太小了 稍微取大一些   10     禁止出现大于等于
    end
    end
end

if iter~=0
    for t=1:si.Horizon
    for s=1:si.NumGl-si.NumGC
        obj4=obj4+penaltyfactor*slack(s,t);
    end
    end
end

for t=1:si.Horizon
    obj5=obj5+0.5*sum(P_P2G(si.P2Glocat_E,t))+si.ELcost*onoff(:,t)+si.stcost*start(:,t);%EL运行成本
%     obj5=obj5+0.5*sum(P_P2G(si.P2Glocat_E,t));
end

for t=1:si.Horizon
    obj6=obj6+0.5*sum(P_de(si.P2Glocat_E,t));%海水淡化运行成本
end
for t=1:si.Horizon
    obj7=obj7+0.5*sum(GMR(si.MRlocat_G,t));%MR维护成本0.2
end
for t=1:si.Horizon
    obj8=obj8+si.Wc*sum(PW_1(:,t)-PW(:,t));
end

for t=1:si.Horizon
    obj9=obj9+si.RBcost*sum(HRB(si.HTlocat_E,t))+si.EBcost*sum(PEB(si.EBlocat,t));%余热回收锅炉和电锅炉成本
end
for t=1:si.Horizon
    obj9a=obj9a+si.CCScost*sum(PCCS(si.CCSlocat,t));%CCS运行成本
end
obj9a=obj9a+0.1*C_STO;%CCS储存成本
% 
obj9b=C_cost;%阶梯式碳交易成本
obj=obj1+obj2+obj3+obj4+obj5+obj6+obj7+obj8+obj9+obj9a+obj9b;
ops=sdpsettings('verbose',1,'solver','gurobi','gurobi.MIPgap',5e-2);
sol=optimize(Constraints,obj,ops) 

ZZ=sdpvar(si.NumGl,si.Horizon);ZZZ=sdpvar(si.NumGl,si.Horizon);
for t=1:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        W=si.TPD(s,4);
        ZZ(s,t)=PHI(s,t)-(Gflow(m,n,t)^2/W^2);%中间变量 检验差值
        ZZZ(s,t)=Prs_square(m,t)-Prs_square(n,t)-(Gflow(m,n,t)^2/W^2);
    end
end

fea_flag = 1;
if sol.problem~=0
    warning(sol.info)
    fea_flag = 0;
else
    mas.slack=double(slack);
    mas.fea_flag=double(fea_flag);
    mas.ZZ=double(ZZ);mas.ZZZ=double(ZZZ);
    mas.onoff=double(onoff);mas.start=double(start);mas.shut=double(shut);
    mas.Gflow=double(Gflow);%支路气流
    mas.PHT_E=double(PHT_E);
    mas.PW_g=double(PW_g);mas.PV_g=double(PV_g);
    mas.PW_cp=double(PW_cp);mas.PV_cp=double(PV_cp);
    mas.P_P2G=double(P_P2G);
    mas.PCCS=double(PCCS);
    mas.GW=double(GW);
    mas.Gflow_1=double(Gflow_1);mas.Gflow_2=double(Gflow_2);mas.Gflow_in1=double(Gflow_in1);mas.Gflow_in2=double(Gflow_in2);mas.Gflow_out1=double(Gflow_out1);mas.Gflow_out2=double(Gflow_out2);
    mas.GHT=double(GHT);mas.GC=double(GC);mas.GP2G=double(GP2G);mas.GHNG=double(GHNG);mas.GHys_in=double(GHys_in);mas.GHys_out=double(GHys_out);mas.pia=double(pia);
    mas.Linepack=double(Linepack);
    mas.uHysch=double(uHysch);mas.uHysdch=double(uHysdch);mas.SHys=double(SHys);
    mas.C_HT=double(C_HT);mas.C_MR=double(C_MR);
    mas.C_q=double(C_q);mas.C_IES=double(C_IES);mas.C_CCS=double(C_CCS);mas.C_a=double(C_a);
    mas.GHMR=double(GHMR);mas.GMR=double(GMR);mas.PMR=double(PMR);
    mas.C_CCS1=double(C_CCS1);mas.C_CCS2=double(C_CCS2);mas.C_cost=double(C_cost);
    mas.v_poor=double(v_poor);mas.V_poor=double(V_poor);mas.v_rich=double(v_rich);mas.V_rich=double(V_rich);
    mas.obj=double(obj);mas.obj1=double(obj1);mas.obj2=double(obj2);mas.obj3=double(obj3);mas.obj4=double(obj4);mas.obj5=double(obj5);mas.obj6=double(obj6);mas.obj7=double(obj7);mas.obj8=double(obj8);mas.obj9=double(obj9);mas.obj9a=double(obj9a);mas.obj9b=double(obj9b);
    mas.PHI=double(PHI);mas.Prs_square=double(Prs_square);
    mas.Pb=double(Pb);mas.Ps=double(Ps);mas.QHT=double(QHT_E);mas.U=double(U);
    mas.PL=double(si.PL);mas.PW=double(PW);
    mas.Ratio=double(Ratio);mas.H_node=double(H_node);
    mas.W1=double(W1);mas.W2=double(W2);
    mas.HRB=double(HRB);mas.PEB=double(PEB);mas.HEB=double(HEB);
%     mas.Hch=double(Hch);mas.Hdch=double(Hdch);mas.HS=double(HS);mas.H=double(H);mas.HE=double(HE);mas.uH=double(uH);
    mas.ES=double(ES);mas.Pdch=double(Pdch);mas.Pch=double(Pch);mas.uEch=double(uEch);mas.uEdch=double(uEdch);
    mas.omiga1=double(omiga1);mas.omiga2=double(omiga2);mas.omiga3=double(omiga3);mas.omiga4=double(omiga4);
    mas.B1=double(B1);mas.B2=double(B2);mas.B3=double(B3);mas.eta_H2=double(eta_H2);mas.P_de=double(P_de);
end
