function [cl] = cluster_dp(dist)
%����ض����dc
percent = 2.0;
N = size(dist,1);
position = round((N*(N-1)/2)*percent/100);
%��ȡdist���������󣬶Խ��߾�Ϊ0
tri_dist = triu(dist,1);
s_dist = sort(tri_dist(tri_dist~=0));
dc = s_dist(position);
%����ضϾ���
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);

%����ֲ��ܶȣ����ø�˹�ˣ�
rho = sum(exp(-(dist./dc).^2),2);
%�����о����е����ֵ
max_dist=max(max(dist));
%�� rho ���������У�ordrho���ԭ��rho_sorted��Ž����rho
[rho_sorted,ordrho]=sort(rho,'descend');

% ���� rho ֵ�������ݵ�
delta(ordrho(1))=-1;
%���ľֲ���ֵ�ĵ㣬nneighΪ0
nneigh(ordrho(1))=0;
% ���� delta �� nneigh ����
for k=2:N
   delta(ordrho(k))=max_dist;
   for p=1:k-1
     if(dist(ordrho(k),ordrho(p))<delta(ordrho(k)))
        delta(ordrho(k))=dist(ordrho(k),ordrho(p));
        nneigh(ordrho(k))=ordrho(p);    % ��¼���оֲ���ֵ��ordrho(k)��ĵ��У���������ĵ�ı��  
     end
   end
end
% ���� rho ֵ������ݵ�� delta ֵ
delta(ordrho(1))=max(delta(:));


%����gammaֵ
for i=1:N
  gamma(i)=rho(i)*delta(i);
end
[gamma_sorted,ordgamma]=sort(gamma,'descend');
plot(gamma_sorted(1:50),'o','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','r');
title ('Decision Graph','FontSize',15.0)

%������ͣ������һ������ѡȡ���������q������ִ�����������
%q = input('q=');
q=26;
% ��ʼ�� cluster ����
NCLUST=0;
%cl(i)=j ��ʾ�� i �����ݵ�����ڵ� j �� cluster
for b=1:N
  cl(b)=-1;
end 
% ͳ�����ݵ㣨���������ģ��ĸ���
for n=1:N
  if gamma(n) > gamma_sorted(q)
     NCLUST=NCLUST+1;
     cl(n)=NCLUST; % �� n �����ݵ��ǵ� NCLUST��ľ�������
     icl(NCLUST)=n;%��ӳ��,�� NCLUST �� cluster ������Ϊ�� n �����ݵ�
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);

%���������ݵ����,����ordrho��nneigh
disp('assignation');
for i=1:N
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end
