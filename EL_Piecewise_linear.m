clear all; close all; clc;
warning off
RTS_Data4;

%% 变量定义
PG=sdpvar(si.NumEn,si.Horizon,'full');%nonNGU C720
dta=sdpvar(si.NumEn,si.Horizon,'full');dta(1,:)=zeros(1,si.Horizon);%1节点是平衡节点 %节点相角 C1440 C1464
P=sdpvar(si.NumEn,si.Horizon,'full');%每个电力节点的净功率注入 C2184
PCCS=sdpvar(si.NumEn,si.Horizon,'full');%碳捕集装置耗电量 C2904
PW_g=sdpvar(si.NumEn,si.Horizon,'full');PV_g=sdpvar(si.NumEn,si.Horizon,'full');%风光上网功率 C3624 C4344
PW_cp=sdpvar(si.NumEn,si.Horizon,'full');PV_cp=sdpvar(si.NumEn,si.Horizon,'full');%风光供给P2G的功率 C5064 C5784
PV=si.PV;PW=si.PW;%风光出力
P_P2G=sdpvar(si.NumEn,si.Horizon,'full');%P2G耗电量 C6504
P_P2G_G=sdpvar(si.NumGn,si.Horizon,'full'); %C6984
PMR=sdpvar(si.NumEn,si.Horizon,'full');%MR设备耗电量C7704

pia=sdpvar(si.NumGn,si.Horizon,'full');%气节点气压 C8184
GW=sdpvar(si.NumGn,si.Horizon,'full');%气源出气量 C8664
GC=sdpvar(si.NumGn,si.Horizon,'full');%压缩机 C9144
Gflow=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路气流 C17808
GHT=sdpvar(si.NumGn,si.Horizon,'full');%HT的耗气量 C18288
GP2G=sdpvar(si.NumGn,si.Horizon,'full');%P2G生产氢气 C18768
GHMR=sdpvar(si.NumGn,si.Horizon,'full');%MR消耗的氢气流量 C19248
GMR=sdpvar(si.NumGn,si.Horizon,'full');%MR产生甲烷流量 C19728
GHNG=sdpvar(si.NumGn,si.Horizon,'full');%P2G向气网输送的氢气 C20208
GHys_in=sdpvar(si.NumGn,si.Horizon,'full');GHys_out=sdpvar(si.NumGn,si.Horizon,'full');%储氢罐的流入、流出流量 C20688 C21168
SHys=sdpvar(si.NumGn,si.Horizon,'full');%储氢罐容量 C21648
uHysch=binvar(si.NumGn,si.Horizon,'full');uHysdch=binvar(si.NumGn,si.Horizon,'full');%表征储氢罐充放的01状态标识 C22128 C22608
GL=si.GL;%气负荷

PHT_E=sdpvar(si.NumEn,si.Horizon,'full');%HT的电功率（电网维度） C23328

C_HT=sdpvar(si.NumGn,si.Horizon,'full');C_G=sdpvar(si.NumEn,si.Horizon,'full');%HT和G的碳排放量，限制条件中补充 C23808 C24528
C_CCS1=sdpvar(si.NumEn,si.Horizon,'full');C_CCS2=sdpvar(si.NumEn,si.Horizon,'full');
v_rich=sdpvar(si.NumEn,si.Horizon,'full');v_poor=sdpvar(si.NumEn,si.Horizon,'full');
V_rich=sdpvar(si.NumEn,si.Horizon,'full');V_poor=sdpvar(si.NumEn,si.Horizon,'full');
C_rich=sdpvar(si.NumEn,si.Horizon,'full');
C_MR=sdpvar(si.NumEn,si.Horizon,'full');

Direction=sdpvar(1,1);%C24529
gamma=sdpvar(1,1);%气网方向指示变量%C24530
Prs_square=sdpvar(si.NumGn,si.Horizon,'full');%压力平方25010
PHI=sdpvar(si.NumGl,si.Horizon,'full');%辅助变量 C25466

Gflow_1=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Gflow_2=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路正向、反向气流
Gflow_in1=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Gflow_out1=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路流入、流出气流
Gflow_in2=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Gflow_out2=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%支路反向流入、流出气流
Linepack=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');Linepack_D=sdpvar(si.NumGn-1,si.NumGn,si.Horizon,'full');%管存

C_cost=sdpvar(1,1);%碳交易成本

onoff=binvar(si.NumPG,si.Horizon,'full');%on-off state of units
shut=binvar(si.NumPG,si.Horizon,'full');%shut-down state of units
start=binvar(si.NumPG,si.Horizon,'full');%start-up state of units
% 
% B1=binvar(1,si.Horizon);B2=binvar(1,si.Horizon);B3=binvar(1,si.Horizon);
% omiga0=sdpvar(1,si.Horizon);omiga1=sdpvar(1,si.Horizon);omiga2=sdpvar(1,si.Horizon);omiga3=sdpvar(1,si.Horizon);
%% 各个设备限制条件
Constraints=[];

%% Electricity grid
Constraints = [Constraints,P==-si.PL];
P(si.PGlocat,:)=P(si.PGlocat,:)+PG(si.PGlocat,:);%火电
P(si.HTlocat_E,:)=P(si.HTlocat_E,:)+PHT_E(si.HTlocat_E,:);%HT
P(si.PWlocat,:)=P(si.PWlocat,:)+PW_g(si.PWlocat,:);%风电
P(si.PVlocat,:)=P(si.PVlocat,:)+PV_g(si.PVlocat,:);%光伏
% P(si.P2Glocat_E,:)=P(si.P2Glocat_E,:)-P_P2G(si.P2Glocat_E,:);%-P2G
P(si.CCSlocat,:)=P(si.CCSlocat,:)-PCCS(si.CCSlocat,:);%-CCS
P(si.MRlocat_E,:)=P(si.MRlocat_E,:)-PMR(si.MRlocat_E,:);%-MR

Constraints = [Constraints,-2*pi<=dta<=2*pi];%节点相角
Constraints = [Constraints,si.Y*dta.*100==P];%节点功率平衡等式   100MVA基准值

for t=1:si.Horizon
    Constraints = [Constraints,PHT_E(si.HTlocat_E,t)==0.1485679853*GHT(si.HTlocat_G,t)];%HT发电功率进行转换（其余行进行清零）0.1485679853
    Constraints = [Constraints,PMR(si.MRlocat_E,t)==0.02*GMR(si.MRlocat_G,t)];%MR耗电
    Constraints = [Constraints,C_HT(si.HTlocat_G,t)==0.1*si.lamuda_CH4*si.aerfa_nv_CH4*si.N_in_CH4*GHT(si.HTlocat_G,t)];%HT碳排量
    Constraints = [Constraints,C_G(si.PGlocat,t)==si.a*PG(si.PGlocat,t)];%火电机组碳排量
end
for t=1:si.Horizon
    Constraints = [Constraints,C_CCS1(si.PGlocat,t)==si.absorb*C_G(si.PGlocat,t)];
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
    Constraints = [Constraints,PCCS(si.HTlocat_E,t)==si.lamuda_CC*si.eta_cc*C_CCS2(si.HTlocat_G,t)];%CCS在HT机组处出力
    Constraints = [Constraints,PCCS(si.PGlocat,t)==si.lamuda_CC*si.eta_cc*C_CCS2(si.PGlocat,t)];%CCS在火电机组处出力
end

PHT_E(setdiff(1:si.NumEn,[si.HTlocat_E]),:)=0;%HT功率无关行清零

% for t=1:si.Horizon
%     Constraints = [Constraints,B1(1,t)+B2(1,t)+B3(1,t)==1];
%     Constraints = [Constraints,omiga0(1,t)+omiga1(1,t)+omiga2(1,t)+omiga3(1,t)==1];
%     Constraints = [Constraints,P_P2G(si.P2Glocat_E,t)==omiga0(1,t)*25+omiga1(1,t)*40+omiga2(1,t)*80+omiga3(1,t)*100];
%     Constraints = [Constraints,GP2G(si.P2Glocat_G,t)==omiga0(1,t)*5000+omiga1(1,t)*8500+omiga2(1,t)*17500+omiga3(1,t)*22090];%P2G的分段线性化
% end

for t=1:si.Horizon
    Constraints = [Constraints,-si.data(:,4)*2<=si.M*dta(:,t)./si.data(:,3)*100<=si.data(:,4)*2];%线路潮流

    Constraints = [Constraints,si.P2Gmin<=P_P2G(:,t)<=si.P2Gmax];
    Constraints = [Constraints,0<=P_P2G_G(:,t)<=si.P2G_Gmax];%P2G最大最小功率限制
    Constraints = [Constraints,P_P2G(:,t)==PW_cp(:,t)+PV_cp(:,t)];%P2G装置由风光供电

    Constraints = [Constraints,si.PMRmin<=PMR(:,t)<=si.PMRmax];%MR设备
    Constraints = [Constraints,0<=PW_g(:,t)<=PW(:,t)];%上网功率
    Constraints = [Constraints,0<=PV_g(:,t)<=PV(:,t)];
    Constraints = [Constraints,0<=PW_cp(:,t)<=PW(:,t)];%给P2G供气功率
    Constraints = [Constraints,0<=PV_cp(:,t)<=PV(:,t)];
    Constraints = [Constraints,PW(:,t)==PW_g(:,t)+PW_cp(:,t)];
    Constraints = [Constraints,PV(:,t)==PV_g(:,t)+PV_cp(:,t)];%风光出力限制
    Constraints = [Constraints,si.CCSmin<=PCCS(:,t)<=si.CCSmax];%碳捕集装置耗能上下限
    Constraints = [Constraints,si.PGmin<=PG(si.PGlocat,t)<=si.PGmax];%火电机组的出力上下限
end
% for t=2:si.Horizon
%     Constraints = [Constraints,-si.ud<=PG(si.PGlocat,t)-PG(si.PGlocat,t-1)<=si.ud];%火电机组爬滑坡速率约束
% end

for t=1:si.Horizon
    if t==1
        Constraints=[Constraints,start(:,t)-shut(:,t)==onoff(:,t)-si.inonoff];
        Constraints=[Constraints,start(:,t)+shut(:,t)<=ones(si.NumPG,1)];
    else
        Constraints=[Constraints,start(:,t)-shut(:,t)==onoff(:,t)-onoff(:,t-1)];
        Constraints=[Constraints,start(:,t)+shut(:,t)<=ones(si.NumPG,1)];
    end
end
for t=2:si.Horizon
    for unit=1:si.NumPG
        indicator=onoff(unit,t)-onoff(unit,t-1);
        range=t:min(si.Horizon,t+si.minup(unit)-1);
        Constraints=[Constraints,onoff(unit,range)>=indicator];
        indicator=onoff(unit,t-1)-onoff(unit,t);
        range=t:min(si.Horizon,t+si.mindown(unit)-1);
        Constraints=[Constraints,onoff(unit,range)<=1-indicator];
    end
end
for t=1:si.Horizon
    Constraints=[Constraints,onoff(:,t).*si.PGmin<=PG(si.PGlocat,t)<=onoff(:,t).*si.PGmax];%thermal unit active power limit
end
for t=2:si.Horizon
    Constraints=[Constraints,(-1).*si.ud<=PG(si.PGlocat,t)-PG(si.PGlocat,t-1)<=si.ud];%thermal unit upward/downward limit
end
PG(setdiff(1:si.NumEn,[si.PGlocat]),:)=0;
%% GAS GRID include direction and linepack

for t=1:si.Horizon
    Constraints = [Constraints,si.Hysmin<=SHys(:,t)<=si.Hysmax];%储氢罐容量上下限
    Constraints = [Constraints,si.GHysmin<=GHys_in(:,t)<=si.GHysmax.*uHysch(:,t)];%储氢罐冲氢流量上下限
    Constraints = [Constraints,si.GHysmin<=GHys_out(:,t)<=si.GHysmax.*uHysdch(:,t)];%储氢罐放氢流量上下限
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

H_node=si.HeatCH4*(1-si.aerfa)+si.HeatH2*si.aerfa;%节点热值

d=binvar(2,si.Horizon);
dd=binvar(2,si.Horizon);
ddd=binvar(2,si.Horizon);%implise 中if elseif

for t=2:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        Constraints = [Constraints,Linepack(m,n,t)-Linepack(m,n,t-1)==Gflow_in1(m,n,t)-Gflow_out1(m,n,t)+Gflow_in2(m,n,t)-Gflow_out2(m,n,t)];%相邻断面管存
        Constraints = [Constraints,Linepack_D(m,n,t)==Linepack(m,n,t)-Linepack(m,n,t-1)];
    end
end
for s=1:si.NumGl-si.NumGC
    m=si.TPD(s,2);
    n=si.TPD(s,3);
    Constraints = [Constraints,Linepack(m,n,1)-si.Linepack0==Gflow_in1(m,n,1)-Gflow_out1(m,n,1)+Gflow_in2(m,n,1)-Gflow_out2(m,n,1)];%相邻断面管存初值
    Constraints = [Constraints,Linepack_D(m,n,1)==Linepack(m,n,1)-si.Linepack0];
end
for t=1:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        Constraints = [Constraints,Linepack(m,n,t)>=si.Linepack0];%管存约束
        Constraints = [Constraints,Linepack_D(m,n,t)>=0];%是否必须大于0
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
%         Constraints = [Constraints,sum(d(:,t))==1];
%         Constraints = [Constraints,implies(d(1,t),[pia(m,t) >= pia(n,t), Direction == 1])];
%         Constraints = [Constraints,implies(d(2,t),[pia(m,t) <= pia(n,t), Direction == 0])];%添加气体流动限制方向
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

%         Constraints = [Constraints,sum(dd(:,t))==1];
%         Constraints = [Constraints,implies(dd(1,t),[pia(m,t) >= pia(n,t), gamma == 1])];
%         Constraints = [Constraints,implies(dd(2,t),[pia(m,t) <= pia(n,t), gamma == -1])];%添加气体流动限制方向

        Constraints = [Constraints,PHI(s,t)>=Gflow(m,n,t)^2/W^2];% A
        Constraints = [Constraints,(1-gamma)/2*(si.pimin(m)^2-si.pimax(n)^2) <= Prs_square(m,t)-Prs_square(n,t)<=(1+gamma)/2*(si.pimax(m)^2-si.pimin(n)^2)];%E
        Constraints = [Constraints,PHI(s,t)>=Prs_square(n,t)-Prs_square(m,t)+(gamma+1)*(si.pimin(m)^2-si.pimax(n)^2)];
        Constraints = [Constraints,PHI(s,t)>=Prs_square(m,t)-Prs_square(n,t)+(gamma-1)*(si.pimax(m)^2-si.pimin(n)^2)];% B
        Constraints = [Constraints,PHI(s,t)<=Prs_square(n,t)-Prs_square(m,t)+(gamma+1)*(si.pimax(m)^2-si.pimin(n)^2)];
        Constraints = [Constraints,PHI(s,t)<=Prs_square(m,t)-Prs_square(n,t)+(gamma-1)*(si.pimin(m)^2-si.pimax(n)^2)];% C  节点压力 数值

        Constraints = [Constraints,Linepack(m,n,t)==0.6*(pia(m,t)+pia(n,t))/2];%管存 管存常数0.6
    end
    for i=1:si.NumGn
        if ismember(i,si.C_in)
            k=find(i==si.C_in);
            Constraints = [Constraints,si.C_ratiomin(k)*pia(si.C_in(k),t)<=pia(si.C_out(k),t)<=si.C_ratiomax(k)*pia(si.C_in(k),t)];%压缩机支路气压  A  (如果加上这个约束的话就相当于默认了gamma==1，既由输入节点流向输出节点)           
            Constraints = [Constraints,Gflow(si.C_in(k),si.C_out(k),t)==(1+si.tao)*GC(k,t)];%压缩机损耗
            Constraints = [Constraints,0<=GC(k,t)<=si.GCmax(k)];%压缩机出力限制
            Constraints = [Constraints,Linepack(si.C_in(k),si.C_out(k),t)==0];%压缩机管道管存

%             Constraints = [Constraints,sum(ddd(:,t))==1];
%             Constraints = [Constraints,implies(ddd(1,t),[pia(si.C_in(k),t) <= pia(si.C_out(k),t), gamma == 1])];
%             Constraints = [Constraints,implies(ddd(2,t),[pia(si.C_in(k),t) >= pia(si.C_out(k),t), gamma == -1])];%添加气体流动限制方向 D
%             Constraints = [Constraints,si.C_ratiomin(k)^2*Prs_square(si.C_in(k),t)+((1-gamma)/2) *(si.pimin(si.C_out(k))^2-si.pimax(si.C_in(k))^2*si.C_ratiomin(k))...
%                                        <= Prs_square(si.C_out(k),t) <= si.C_ratiomax(k)^2*Prs_square(si.C_in(k),t)+((1-gamma)/2) *(si.pimax(si.C_out(k))^2-si.pimin(si.C_in(k))^2*si.C_ratiomax(k))];%压缩机支路气压 B
%             Constraints = [Constraints,si.C_ratiomin(k)^2*Prs_square(si.C_out(k),t)+((1+gamma)/2) *(si.pimin(si.C_in(k))^2-si.pimax(si.C_out(k))^2*si.C_ratiomin(k))...
%                                        <= Prs_square(si.C_in(k),t) <= si.C_ratiomax(k)^2*Prs_square(si.C_out(k),t)+((1+gamma)/2) *(si.pimax(si.C_in(k))^2-si.pimin(si.C_out(k))^2*si.C_ratiomax(k))];%压缩机支路气压 C
        end
    end
    GHT(setdiff(1:si.NumGn,[si.HTlocat_G]),:)=0;%非HT和P2G节点的耗气量、产气量置为零
    GHNG(setdiff(1:si.NumGn,[si.P2Glocat_G]),:)=0;
    GMR(setdiff(1:si.NumGn,[si.MRlocat_G]),:)=0;
    for i=1:si.NumGn
        if ismember(i,si.P2Glocat_G)
            kk=find(i==si.P2Glocat_G);

            Constraints = [Constraints,GP2G(i,t)==(P_P2G(si.P2Glocat_E(kk),t)*si.eta_H2*3.6*1000*0.5)/si.HeatH2];%P2G电-气转换 10 4
            Constraints = [Constraints,GMR(i,t)==(4*si.eta_MR*GHMR(si.MRlocat_G(kk),t )*si.HeatH2)/si.HeatCH4];%MR产甲烷
            Constraints = [Constraints,GHNG(i,t)==GP2G(si.P2Glocat_G(kk),t)-GHMR(si.P2Glocat_G(kk),t)-GHys_in(si.P2Glocat_G(kk),t)];
        end
    end

    temp=sdpvar(si.NumGn,1);
    for i=1:si.NumGn
        temp(i)=GW(i,t)*si.HeatCH4+GMR(i,t)*si.HeatCH4+GHNG(i,t)*si.HeatH2+GHys_out(i,t)*si.HeatH2-si.GL(i,t)*H_node-GHT(i,t)*H_node;%天然气节点净注入气流   维数
        for s=1:si.NumGl-si.NumGC
            m=si.TPD(s,2);
            n=si.TPD(s,3);
            Constraints = [Constraints,sum(d(:,t))==1];
            Constraints = [Constraints,implies(d(1,t),[pia(m,t) >= pia(n,t), Direction == 1])];
            Constraints = [Constraints,implies(d(2,t),[pia(m,t) <= pia(n,t), Direction == 0])];%添加气体流动限制方向
%             Direction = 1;
            if i==m
                temp(i)=temp(i)-Direction*Gflow(m,n,t)*H_node+(1-Direction)*Gflow(m,n,t)*H_node-Linepack(m,n,t);
%                 temp(i)=temp(i)-Gflow(m,n,t);
            elseif i==n
                temp(i)=temp(i)+Direction*Gflow(m,n,t)*H_node-(1-Direction)*Gflow(m,n,t)*H_node+Linepack(m,n,t);
%                 temp(i)=temp(i)+Gflow(m,n,t);
            end
        end
        if ismember(i,si.C_in)
            k=find(i==si.C_in);
            temp(i)=temp(i)-Gflow(si.C_in(k),si.C_out(k),t)*H_node;
        end
        if ismember(i,si.C_out)
            kk=find(i==si.C_out);
            temp(i)=temp(i)+Gflow(si.C_in(kk),si.C_out(kk),t)*H_node;
        end
        Constraints = [Constraints,temp(i)==0];%天然气节点气流平衡
    end
end

for t=1:si.Horizon
    Constraints = [Constraints,C_MR(si.MRlocat_E,t)==2*(GMR(si.MRlocat_G,t)/(1.397*1e-3))*44*1e-6];%MR消耗的CO2
    Constraints = [Constraints,0<=sum(C_MR(si.MRlocat_E,t))<=sum(si.eta_cc*C_CCS2(si.CCSlocat,t))];
    Constraints = [Constraints,sum(GHNG(si.P2Glocat_G,t))+sum(GHys_out(si.P2Glocat_G,t))==0.1*(sum(GW(si.GWlocat,t))+sum(GMR(si.MRlocat_G,t)))];
end

C_STO=0;
for t=1:si.Horizon
    C_STO=C_STO+sum(si.eta_cc*C_CCS2(si.CCSlocat,t))-sum(C_MR(si.MRlocat_E,t));
end
%% Carbon trade
C_q=0;C_IES=0;C_a=0;C_CCS=0;
for t=1:si.Horizon
    C_q=C_q+si.daita_e*(sum(PG(si.PGlocat,t))+sum(PHT_E(si.HTlocat_E,t)));%碳配额
end

for t=1:si.Horizon
    C_IES=C_IES+(sum(C_G(si.PGlocat,t))+sum(C_HT(si.HTlocat_G,t))+sum((1-si.aerfa)*2.165*1e-3*si.GL(:,t)));%全部碳排
end

for t=1:si.Horizon
    C_CCS=C_CCS+si.eta_cc*(sum(C_G(si.PGlocat,t))+sum(C_HT(si.HTlocat_G,t)));%碳捕集量
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
obj1=0;obj2=0;obj3=0;obj4=0;obj5=0;obj6=0;obj7=0;
for t=1:si.Horizon
    obj1=obj1+PG(si.PGlocat,t)'*diag(si.NGa)*PG(si.PGlocat,t)+si.NGb'*PG(si.PGlocat,t)+sum(si.NGc)+si.c*onoff(:,t)+si.stcost*start(:,t)+si.sdcost*shut(:,t);%PG燃煤成本
%       obj1=obj1+PG(si.PGlocat,t)'*diag(si.NGa)*PG(si.PGlocat,t)+si.NGb'*PG(si.PGlocat,t)+sum(si.NGc);
end
for t=1:si.Horizon
    obj2=obj2+si.fp*sum(GW(si.GWlocat,t));%气源买气成本
end
for t=1:si.Horizon
    obj3=obj3+si.HTcost*sum(PHT_E(si.HTlocat_E,t));%HT燃料成本
end
for t=1:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        obj4=obj4+0.1*PHI(s,t);%weymouth等式松弛惩罚项   0.1取得太小了 稍微取大一些   10     禁止出现大于等于
    end
end
for t=1:si.Horizon
    obj5=obj5+si.CCScost*sum(PCCS(si.CCSlocat,t));%CCS运行成本
end
obj5=obj5+0.1*C_STO;%CCS储存成本

obj6=C_cost;%阶梯式碳交易成本
for t=1:si.Horizon
    obj7=obj7+si.MRcost*sum(GMR(si.MRlocat_G,t));%MR维护成本
end
obj=obj1+obj2+obj3+obj4+obj5+obj6+obj7;
% ops = sdpsettings('verbose',1,'solver','ipopt', 'debug',1,'usex0',0);
ops=sdpsettings('verbose',1,'solver','gurobi','gurobi.MIPgap',5e-3);
sol=optimize(Constraints,obj,ops)   
% for t=1:si.Horizon
%     PFlow(:,t)=si.M*double(dta(:,t))./si.data(:,3)*100;%支路潮流
% end
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


ZZ=double(ZZ);ZZZ=double(ZZZ);
onoff=double(onoff);start=double(start);shut=double(shut);
temp=double(temp);
Gflow=double(Gflow);%支路气流
PG=double(PG);PHT_E=double(PHT_E);
PW_g=double(PW_g);PV_g=double(PV_g);
PW_cp=double(PW_cp);PV_cp=double(PV_cp);
P_P2G=double(P_P2G);
PCCS=double(PCCS);
GW=double(GW);
Gflow_1=double(Gflow_1);Gflow_2=double(Gflow_2);Gflow_in1=double(Gflow_in1);Gflow_in2=double(Gflow_in2);Gflow_out1=double(Gflow_out1);Gflow_out2=double(Gflow_out2);
GHT=double(GHT);GC=double(GC);GP2G=double(GP2G);GHNG=double(GHNG);GHys_in=double(GHys_in);GHys_out=double(GHys_out);pia=double(pia);
Linepack=double(Linepack);Linepack_D=double(Linepack_D);
uHysch=double(uHysch);uHysdch=double(uHysdch);SHys=double(SHys);
C_HT=double(C_HT);C_G=double(C_G);C_MR=double(C_MR);
C_q=double(C_q);C_IES=double(C_IES);C_CCS=double(C_CCS);C_a=double(C_a);
GHMR=double(GHMR);GMR=double(GMR);PMR=double(PMR);Direction=double(Direction);
obj=double(obj);obj1=double(obj1);obj2=double(obj2);obj3=double(obj3);obj4=double(obj4);obj5=double(obj5);obj6=double(obj6);obj7=double(obj7);
PHI=double(PHI);Prs_square=double(Prs_square);
C_CCS1=double(C_CCS1);C_CCS2=double(C_CCS2);C_cost=double(C_cost);
v_poor=double(v_poor);V_poor=double(V_poor);v_rich=double(v_rich);V_rich=double(V_rich);
% omiga0=double(omiga0);omiga1=double(omiga1);omiga2=double(omiga2);omiga3=double(omiga3);
%% carbon emission and gas linepack
C_emission=sdpvar(1,si.Horizon);
for t=1:si.Horizon
C_emission(1,t)=sum(C_HT(si.HTlocat_G,t))+sum(C_G(si.PGlocat,t))+sum((1-si.aerfa)*2.165*1e-3*si.GL(:,t));
end
C_emission=double(C_emission);
C_storage=sdpvar(1,si.Horizon);
for t=1:si.Horizon
    C_storage(1,t)=sum(si.eta_cc*C_CCS2(si.CCSlocat,t))-sum(C_MR(si.MRlocat_E,t));
end
C_storage=double(C_storage);
C_capture=sdpvar(1,si.Horizon);
for t=1:si.Horizon
    C_capture(1,t)=si.eta_cc*(sum(C_G(si.PGlocat,t))+sum(C_HT(si.HTlocat_G,t)));
end
C_capture=double(C_capture);

Gaslinepack=zeros(si.NumGl-si.NumGC,si.Horizon);
for t=1:si.Horizon
    for s=1:si.NumGl-si.NumGC
        m=si.TPD(s,2);
        n=si.TPD(s,3);
        Gaslinepack(s,t)=Linepack(m,n,t);
    end
end
Gaslinepack=double(Gaslinepack);
%% *************重点
% [model,recoverymodel]=export(Constraints,obj,ops);
% iis=gurobi_iis(model);
% gurobi_write(model,'TestModel.lp');
% model.computeIIS();
%% Plot
figure ('NumberTitle','off','Name','Power');
y=[sum(PG(si.PGlocat,:));sum(PHT_E(si.HTlocat_E,:));sum(si.PW);sum(si.PV);-sum(P_P2G);-sum(PCCS);-sum(PMR)];
y(y<1&y>-1)=0;
yy=[];
for t=1:si.Horizon
    yy=[yy;y(:,t)'];
end
t1 = yy;
t2 = yy;
t1(t1<0) = nan;
t2(t2>0) = nan;
c=bar(t1(:,[1:4]),'stacked');
set(c(1),'FaceColor',[0.3 0.5470 0.8410],'EdgeColor','none');
set(c(2),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none');
set(c(3),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none');
set(c(4),'FaceColor',[0.6290 0.740 0.8250],'EdgeColor','none');
hold on
c=bar(t2(:,[5:7]),'stacked');
set(c(1),'FaceColor',[0.6290 0.10 0.4250],'EdgeColor','none');
set(c(2),'FaceColor',[0.4290 0.240 0.2250],'EdgeColor','none');
set(c(3),'FaceColor',[0.1290 0.440 0.6250],'EdgeColor','none');
% legend('PG','HT','PW','PV','P2G','PCCS','PMR');
% legend('boxoff')
xlabel('Time (h)');
ylabel('Power (MW)');
xticks([0 5 10 15 20 25]);
set(gca,'FontSize',20,'Fontname','times new Roman');
grid on


figure ('NumberTitle','off','Name','Gas');
y=[sum(GW);sum(GHMR);sum(GHNG(si.P2Glocat_G,:));sum(GHys_out);-sum(GHT);-sum(GHys_in)];
y(y<1&y>-1)=0;
yy=[];
for t=1:si.Horizon
    yy=[yy;y(:,t)'];
end
t1 = yy;
t2 = yy;
t1(t1<0) = nan;
t2(t2>0) = nan;
c=bar(t1(:,[1:4]),'stacked');
set(c(1),'FaceColor',[0 0.2470 0.5410],'EdgeColor','none');
set(c(2),'FaceColor',[0.9 0.9 0],'EdgeColor','none');
set(c(3),'FaceColor',[0.0000 0.5250 0.0980],'EdgeColor','none');
set(c(4),'FaceColor',[148,0,211]/255,'EdgeColor','none');
% set(c(4),'FaceColor',[0.3000 0.8200 0.0980],'EdgeColor','none');
hold on
c=bar(t2(:,[5:6]),'stacked');
% set(c(1),'FaceColor',[0.6000 0.0000 0.20],'EdgeColor','none');
set(c(1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none');
set(c(2),'FaceColor',[0.000 0.9 0.9],'EdgeColor','none');
% set(c(4),'FaceColor',[0.4290 0.5940 0.1250],'EdgeColor','none');
% legend('GW','GHMR','GHNG','GHysout','GL','GHT','GHysin','GCOM');
% legend('boxoff')
xlabel('Time (h)');
ylabel('Gas (m^3)');
xticks([0 5 10 15 20 25]);
set(gca,'FontSize',20,'Fontname','times new Roman');
grid on



figure ('NumberTitle','off','Name','carbon2');
fig = figure;                               %使用自定义变量fig拿到figure的handle
left_color = [0    0    0];  %百度到的MATLAB默认的铜橙色
right_color = [0  0  1];      %百度到的MATLAB默认的浅蓝色
set(fig,'defaultAxesColorOrder',[left_color; right_color]);   %设置左右纵轴的颜色
x=1:1:si.Horizon;%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
y=[];yy=[];yyy=[];yyyy=[];
for t=1:si.Horizon
y=[y,sum(PG(si.PGlocat,t))/(sum(PG(si.PGlocat,t))+sum(PHT_E(si.HTlocat_E,t))+sum(PV(si.PVlocat,t))+sum(PW(si.PWlocat,t)))];
end
yyyy=[C_emission;C_capture;C_storage;sum(C_MR(si.MRlocat_E,:))];
yyyy(yyyy<1&yyyy>-1)=0;
aaa=[];
for t=1:si.Horizon
    aaa=[aaa;yyyy(:,t)'];
end
t1 = aaa;
t2 = aaa;
t1(t1<0) = nan;
c=bar(t1(:,[1:4]),'stacked');
set(c(1),'FaceColor',[0.8500 0.3250 0.3980],'EdgeColor','none');
set(c(2),'FaceColor',[0.3 0.7470 0.8410],'EdgeColor','none');
set(c(3),'FaceColor',[0.9290 0.6940 0.3250],'EdgeColor','none');
set(c(4),'FaceColor',[0.3290 0.540 0.2250],'EdgeColor','none');
hold on
xlabel('Time (h)')  %x轴坐标描述
ylabel('Carbon (t)') %y轴坐标描述

yyaxis right;
plot(x,y,'-p','color','b','MarkerSize',8,'LineWidth',3);
hold on   %x轴坐标描述
ylabel('Proportion','Color','b') %y轴坐标描述

xticks([0 5 10 15 20 25]);
% legend('Carbon emission','Carbon capture','Carbon storage','MR consume carbon','Proportion of coal-fired units');   %右上角标注
set(gca,'FontSize',20,'Fontname','times new Roman');
grid on

figure ('NumberTitle','off','Name','Carbon3');
fig = figure;                               %使用自定义变量fig拿到figure的handle
left_color = [0    0    0];  %百度到的MATLAB默认的铜橙色
right_color = [0    0    1];      %百度到的MATLAB默认的浅蓝色
x=1:1:si.Horizon;%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
y=[];yy=[];yyy=[];yyyy=[];yyyyy=[];oo=zeros(1,24);
for t=1:si.Horizon
y=[y,sum(v_rich(si.CCSlocat,t))];
end
yyyy=[sum(V_rich(si.CCSlocat,:));sum(V_poor(si.CCSlocat,:))];
yyyy(yyyy<1&yyyy>-1)=0;
% figure ('NumberTitle','off','Name','Carbon3');
% fig = figure;                               %使用自定义变量fig拿到figure的handle
% left_color = [0    0    0];  %百度到的MATLAB默认的铜橙色
% right_color = [0    0    0];      %百度到的MATLAB默认的浅蓝色
x=1:1:si.Horizon;%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
aaa=[];
for t=1:si.Horizon
    aaa=[aaa;yyyy(:,t)'];
end
t1 = aaa;
t2 = aaa;
t1(t1<0) = nan;
c=bar(t1(:,[1:2]));
% yyaxis left;
set(c(1),'FaceColor',[0.3290 0.6940 0.1250],'EdgeColor','none');
set(c(2),'FaceColor',[0.3290 0.540 0.8250],'EdgeColor','none');
hold on
xlabel('Time (h)')  %x轴坐标描述
ylabel('Volume (m^3)') %y轴坐标描述
% yyaxis right;
plot(x,y,'-p','color','r','MarkerSize',8,'LineWidth',5);
hold on   %x轴坐标描述
xticks([0 5 10 15 20 25]);
% legend('Volume of rich liquid tank','Volume of poor liquid tank','Flow of rich liquid tank');   %右上角标注
set(gca,'FontSize',32,'Fontname','times new Roman');
grid on
