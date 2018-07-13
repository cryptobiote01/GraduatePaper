function topweijiaqian(nelx,nely,volfrac,penal,rmin);  
  
nelx=40;        % x轴方向上单元个数  
nely=40;        % y轴方向上单元个数  
volfrac=0.6;     %体积比  
penal=3;        %材料插值的惩罚因子  
rmin=2;         %敏度过滤半径  
a = 0.003;         %多目标函数权重
ks=0.001;

  
% INITIALIZE  
x(1:nely,1:nelx) = volfrac;  % x是设计变量(单元伪密度)  
loop = 0;                   %存放迭代次数的变量  
change = 1.;                    %每次迭代，目标函数（柔度）的改变值，用来判断何时收敛  


num=0
min=[];
  
% START ITERATION  
while change > 0.01 %当两次连续目标函数迭代的差<=0.01时，迭代结束  
  loop = loop + 1;  
  xold = x; %把前一次的设计变量付给xold  
  
% FE-ANALYSIS  
  [U]=FE(nelx,nely,x,penal); %有限元分析，得到位移矢量U  
  [V]=FE1(nelx,nely,x,penal);
  
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS  
  [KE] = lk;     %单位刚度矩阵  
  c = 0.;       %用来存放目标函数的变量，这里刚度最大，柔度最小  
  
  for ely = 1:nely  
    for elx = 1:nelx  
      n1 = (nely+1)*(elx-1)+ely;  %左上角的单元节点  
      n2 = (nely+1)* elx   +ely;  %右上角的单元节点  
  
      %所示单元的自由度，左上，右上，右下，左下  
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);  
      Ve = V([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);  
  
%       c = c + x(ely,elx)^penal*(-Ve'*KE*Ue)./Ue'*KE*Ue; %计算目标函数的值（柔度）  
%       dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*(-Ve'*KE*Ue)./Ue'*KE*Ue; %目标函数的灵敏度  
      
%       c = c + x(ely,elx)^penal*(-a*(Ve'*KE*Ue)+(1-a)*Ue'*KE*Ue); %计算目标函数的值（柔度）  
%       dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*(-a*(Ve'*KE*Ue)+(1-a)*Ue'*KE*Ue); %目标函数的灵敏度  

        c = c + x(ely,elx)^penal*ks*(Ve'*KE*Ue)^2/(2*Ue'*KE*Ue);
        dc(ely,elx) =-penal*x(ely,elx)^(penal-1)*ks*(Ve'*KE*Ue)^2/(2*Ue'*KE*Ue); %目标函数的灵敏度  
      %dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; %目标函数的灵敏度  
    end  
  end  
  
% FILTERING OF SENSITIVITIES  
  [dc]   = check(nelx,nely,rmin,x,dc); %灵敏度过滤，为了边界光顺一点   
  
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD  
  [x]    = OC(nelx,nely,x,volfrac,dc);   
  
% PRINT RESULTS 屏幕上显示迭代信息  
  change = max(max(abs(x-xold))); %计算目标函数的改变量  
%   disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...  
%        ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...  
%         ' ch.: ' sprintf('%6.3f',change )])  
  
% PLOT DENSITIES 优化结果的图形显示   
%x=abs(x)
%     min=[min c]
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);  
% %   num=num+1;
% %   name = ['save' num2str(num) '.jpg'];
% %   saveas(gcf,name);
%     if loop==150
%         minn=sort(min,'descend');
%         plot(minn,'k','LineWidth',2);
%         xlabel('迭代次数');
%         ylabel('目标函数值');
%         saveas(gcf,'min.jpg');
%         break;
%     end
end   
  
  
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [xnew]=OC(nelx,nely,x,volfrac,dc)    
  
l1 = 0; l2 = 100000; %用于体积约束的拉格朗日乘子  
move = 0.2;  
  
while (l2-l1 > 1e-4)  
  lmid = 0.5*(l2+l1);  
  
  %即论文公式的综合  
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));  
  
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;%二乘法减半  
    l1 = lmid;  
  else  
    l2 = lmid;  
  end  
end  
  
  
%%%%%%%%%% MESH-INDEPENDENCY FILTER 敏度过滤技术子程序%%%%%%%%%%%%%%%%%%%  
function [dcn]=check(nelx,nely,rmin,x,dc)  
  
dcn=zeros(nely,nelx);  
  
for i = 1:nelx  
  for j = 1:nely  
sum=0.0;   
  
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)  
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)  
        fac = rmin-sqrt((i-k)^2+(j-l)^2);  
        sum = sum + max(0,fac);  
  
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);  
      end  
end  
  
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);  
  end  
end  
  
  
%%%%%%%%%% FE-ANALYSIS 有限元求解子程序%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [U]=FE(nelx,nely,x,penal) %自定义函数，最后返回[U]  
  
[KE] = lk; %单元刚度矩阵  
  
% sparse 把一个全矩阵转化为一个稀疏矩阵，只存储每一个非零元素的三个值：元素值，元素的行号和列号  
%总体刚度矩阵的稀疏矩阵  
% *2是因为x，y都有一个数  
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));  
%力矩阵的稀疏矩阵  
F = sparse(2*(nely+1)*(nelx+1),1);   
U = zeros(2*(nely+1)*(nelx+1),1); %零矩阵  
  
for elx = 1:nelx  
  for ely = 1:nely  
    %一列列的排序  
    n1 = (nely+1)*(elx-1)+ely; %左上  
n2 = (nely+1)* elx   +ely;  %右上  
  
% 左上，右上，右下，左下 自由度  
% 一个点有两个，所以要*2。第一个从1开始，所以*2之后要-1。  
edof = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];  
  
%将单元刚度矩阵组装成总的刚度矩阵  
K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;    
end  
end  
  
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)  
F(819,1) = 100; % 应用了一个在左上角的垂直单元力。  
%按着图上来的，最左边和右下角已经固定  
fixeddofs   = union([1:2*(nely+1)],[820]); %固定结点  
alldofs     = [1:2*(nely+1)*(nelx+1)];                      %所有结点  
  
% setdiff 因无约束自由度与固定自由度的不同来找到无约束自由度  
freedofs    = setdiff(alldofs,fixeddofs); %不受约束的自由度  
  
% SOLVING  
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);        
U(fixeddofs,:)= 0; % 矩阵A的第r行：A（r，：）   


function [U]=FE1(nelx,nely,x,penal) %自定义函数，最后返回[U]  
  
[KE] = lk; %单元刚度矩阵  
  
% sparse 把一个全矩阵转化为一个稀疏矩阵，只存储每一个非零元素的三个值：元素值，元素的行号和列号  
%总体刚度矩阵的稀疏矩阵  
% *2是因为x，y都有一个数  
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));  
%力矩阵的稀疏矩阵  
F = sparse(2*(nely+1)*(nelx+1),1);   
U = zeros(2*(nely+1)*(nelx+1),1); %零矩阵  
  
for elx = 1:nelx  
  for ely = 1:nely  
    %一列列的排序  
    n1 = (nely+1)*(elx-1)+ely; %左上  
n2 = (nely+1)* elx   +ely;  %右上  
  
% 左上，右上，右下，左下 自由度  
% 一个点有两个，所以要*2。第一个从1开始，所以*2之后要-1。  
edof = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];  
  
%将单元刚度矩阵组装成总的刚度矩阵  
K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;    
end  
end  
  
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)  
F(1934,1) = -3; % 应用了一个在左上角的垂直单元力。  
%按着图上来的，最左边和右下角已经固定  
fixeddofs   = union([1:2*(nely+1)],[820]); %固定结点  
alldofs     = [1:2*(nely+1)*(nelx+1)];                      %所有结点  
  
% setdiff 因无约束自由度与固定自由度的不同来找到无约束自由度  
freedofs    = setdiff(alldofs,fixeddofs); %不受约束的自由度  
  
% SOLVING  
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);        
U(fixeddofs,:)= 0; % 矩阵A的第r行：A（r，：）   



  
%%%%%%%%%% ELEMENT STIFFNESS MATRIX 单元刚度矩阵的子程序%%%%%%%%%%%%%%%%%%%%  
function [KE]=lk  
  
E =5.;   
nu = 0.3;  
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...   
-1/4+nu/12  -1/8-nu/8  nu/6       1/8-3*nu/8];  
  
                %u1,v1,  u2,v2,   u3,v3,   u4,v4  
KE = E/(1-nu^2)*[    k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)  
                     k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)  
                     k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)  
                     k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)  
                     k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)  
                     k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)  
                     k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)  
                     k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];  