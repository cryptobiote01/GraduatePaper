function topweijiaqian(nelx,nely,volfrac,penal,rmin);  
  
nelx=40;        % x�᷽���ϵ�Ԫ����  
nely=40;        % y�᷽���ϵ�Ԫ����  
volfrac=0.6;     %�����  
penal=3;        %���ϲ�ֵ�ĳͷ�����  
rmin=2;         %���ȹ��˰뾶  
a = 0.003;         %��Ŀ�꺯��Ȩ��
ks=0.001;

  
% INITIALIZE  
x(1:nely,1:nelx) = volfrac;  % x����Ʊ���(��Ԫα�ܶ�)  
loop = 0;                   %��ŵ��������ı���  
change = 1.;                    %ÿ�ε�����Ŀ�꺯������ȣ��ĸı�ֵ�������жϺ�ʱ����  


num=0
min=[];
  
% START ITERATION  
while change > 0.01 %����������Ŀ�꺯�������Ĳ�<=0.01ʱ����������  
  loop = loop + 1;  
  xold = x; %��ǰһ�ε���Ʊ�������xold  
  
% FE-ANALYSIS  
  [U]=FE(nelx,nely,x,penal); %����Ԫ�������õ�λ��ʸ��U  
  [V]=FE1(nelx,nely,x,penal);
  
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS  
  [KE] = lk;     %��λ�նȾ���  
  c = 0.;       %�������Ŀ�꺯���ı���������ն���������С  
  
  for ely = 1:nely  
    for elx = 1:nelx  
      n1 = (nely+1)*(elx-1)+ely;  %���Ͻǵĵ�Ԫ�ڵ�  
      n2 = (nely+1)* elx   +ely;  %���Ͻǵĵ�Ԫ�ڵ�  
  
      %��ʾ��Ԫ�����ɶȣ����ϣ����ϣ����£�����  
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);  
      Ve = V([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);  
  
%       c = c + x(ely,elx)^penal*(-Ve'*KE*Ue)./Ue'*KE*Ue; %����Ŀ�꺯����ֵ����ȣ�  
%       dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*(-Ve'*KE*Ue)./Ue'*KE*Ue; %Ŀ�꺯����������  
      
%       c = c + x(ely,elx)^penal*(-a*(Ve'*KE*Ue)+(1-a)*Ue'*KE*Ue); %����Ŀ�꺯����ֵ����ȣ�  
%       dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*(-a*(Ve'*KE*Ue)+(1-a)*Ue'*KE*Ue); %Ŀ�꺯����������  

        c = c + x(ely,elx)^penal*ks*(Ve'*KE*Ue)^2/(2*Ue'*KE*Ue);
        dc(ely,elx) =-penal*x(ely,elx)^(penal-1)*ks*(Ve'*KE*Ue)^2/(2*Ue'*KE*Ue); %Ŀ�꺯����������  
      %dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; %Ŀ�꺯����������  
    end  
  end  
  
% FILTERING OF SENSITIVITIES  
  [dc]   = check(nelx,nely,rmin,x,dc); %�����ȹ��ˣ�Ϊ�˱߽��˳һ��   
  
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD  
  [x]    = OC(nelx,nely,x,volfrac,dc);   
  
% PRINT RESULTS ��Ļ����ʾ������Ϣ  
  change = max(max(abs(x-xold))); %����Ŀ�꺯���ĸı���  
%   disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...  
%        ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...  
%         ' ch.: ' sprintf('%6.3f',change )])  
  
% PLOT DENSITIES �Ż������ͼ����ʾ   
%x=abs(x)
%     min=[min c]
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);  
% %   num=num+1;
% %   name = ['save' num2str(num) '.jpg'];
% %   saveas(gcf,name);
%     if loop==150
%         minn=sort(min,'descend');
%         plot(minn,'k','LineWidth',2);
%         xlabel('��������');
%         ylabel('Ŀ�꺯��ֵ');
%         saveas(gcf,'min.jpg');
%         break;
%     end
end   
  
  
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [xnew]=OC(nelx,nely,x,volfrac,dc)    
  
l1 = 0; l2 = 100000; %�������Լ�����������ճ���  
move = 0.2;  
  
while (l2-l1 > 1e-4)  
  lmid = 0.5*(l2+l1);  
  
  %�����Ĺ�ʽ���ۺ�  
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));  
  
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;%���˷�����  
    l1 = lmid;  
  else  
    l2 = lmid;  
  end  
end  
  
  
%%%%%%%%%% MESH-INDEPENDENCY FILTER ���ȹ��˼����ӳ���%%%%%%%%%%%%%%%%%%%  
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
  
  
%%%%%%%%%% FE-ANALYSIS ����Ԫ����ӳ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [U]=FE(nelx,nely,x,penal) %�Զ��庯������󷵻�[U]  
  
[KE] = lk; %��Ԫ�նȾ���  
  
% sparse ��һ��ȫ����ת��Ϊһ��ϡ�����ֻ�洢ÿһ������Ԫ�ص�����ֵ��Ԫ��ֵ��Ԫ�ص��кź��к�  
%����նȾ����ϡ�����  
% *2����Ϊx��y����һ����  
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));  
%�������ϡ�����  
F = sparse(2*(nely+1)*(nelx+1),1);   
U = zeros(2*(nely+1)*(nelx+1),1); %�����  
  
for elx = 1:nelx  
  for ely = 1:nely  
    %һ���е�����  
    n1 = (nely+1)*(elx-1)+ely; %����  
n2 = (nely+1)* elx   +ely;  %����  
  
% ���ϣ����ϣ����£����� ���ɶ�  
% һ����������������Ҫ*2����һ����1��ʼ������*2֮��Ҫ-1��  
edof = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];  
  
%����Ԫ�նȾ�����װ���ܵĸնȾ���  
K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;    
end  
end  
  
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)  
F(819,1) = 100; % Ӧ����һ�������ϽǵĴ�ֱ��Ԫ����  
%����ͼ�����ģ�����ߺ����½��Ѿ��̶�  
fixeddofs   = union([1:2*(nely+1)],[820]); %�̶����  
alldofs     = [1:2*(nely+1)*(nelx+1)];                      %���н��  
  
% setdiff ����Լ�����ɶ���̶����ɶȵĲ�ͬ���ҵ���Լ�����ɶ�  
freedofs    = setdiff(alldofs,fixeddofs); %����Լ�������ɶ�  
  
% SOLVING  
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);        
U(fixeddofs,:)= 0; % ����A�ĵ�r�У�A��r������   


function [U]=FE1(nelx,nely,x,penal) %�Զ��庯������󷵻�[U]  
  
[KE] = lk; %��Ԫ�նȾ���  
  
% sparse ��һ��ȫ����ת��Ϊһ��ϡ�����ֻ�洢ÿһ������Ԫ�ص�����ֵ��Ԫ��ֵ��Ԫ�ص��кź��к�  
%����նȾ����ϡ�����  
% *2����Ϊx��y����һ����  
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));  
%�������ϡ�����  
F = sparse(2*(nely+1)*(nelx+1),1);   
U = zeros(2*(nely+1)*(nelx+1),1); %�����  
  
for elx = 1:nelx  
  for ely = 1:nely  
    %һ���е�����  
    n1 = (nely+1)*(elx-1)+ely; %����  
n2 = (nely+1)* elx   +ely;  %����  
  
% ���ϣ����ϣ����£����� ���ɶ�  
% һ����������������Ҫ*2����һ����1��ʼ������*2֮��Ҫ-1��  
edof = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];  
  
%����Ԫ�նȾ�����װ���ܵĸնȾ���  
K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;    
end  
end  
  
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)  
F(1934,1) = -3; % Ӧ����һ�������ϽǵĴ�ֱ��Ԫ����  
%����ͼ�����ģ�����ߺ����½��Ѿ��̶�  
fixeddofs   = union([1:2*(nely+1)],[820]); %�̶����  
alldofs     = [1:2*(nely+1)*(nelx+1)];                      %���н��  
  
% setdiff ����Լ�����ɶ���̶����ɶȵĲ�ͬ���ҵ���Լ�����ɶ�  
freedofs    = setdiff(alldofs,fixeddofs); %����Լ�������ɶ�  
  
% SOLVING  
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);        
U(fixeddofs,:)= 0; % ����A�ĵ�r�У�A��r������   



  
%%%%%%%%%% ELEMENT STIFFNESS MATRIX ��Ԫ�նȾ�����ӳ���%%%%%%%%%%%%%%%%%%%%  
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