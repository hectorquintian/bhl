function w=bhl(data,a,b,lRate,iters,m,id,classes)
close all
[fil,col]=size(data);
colors={'b.';'g.';'rd';'c.';'m.';'y.';'k.';'bx';'go';'r.'}; %Lista de colores con los que se ploteará cada clase
tnum=id;
%%
%%%%PREPOSESS%%%%%
% varD=var(data);
% meanD=mean(data);
%  for i=1:size(data,2)
%      data(:,i)=data(:,i)-meanD(i);
%      data(:,i)=data(:,i)/varD(i);
%  end
% data=data/(max(max(data)));
%data=normalizar_datos(data,0,1);
% minD=min(min(data));
% data=data-minD;
% maxD=max(max(data));
% data=data/maxD;
figure
hist(data)

ex=0:0.0001:1;
ey=0:0.0001:1;
e=[ex;ey];
%LEARNING RULE BETA-SIM
figure
title(strcat('BETA-SIM learning rule'))
for i=1:size(a,2)
%    beta=(factorial(a(i)-1)*factorial(b(i)-1)/factorial(a(i)+b(i)-1));
    for j=1:size(e,2)
        err2=e(2,j);

        err2=abs(err2);
        d21=(err2^(a-2)) * ((1-err2)^(b-2));
        d22=((-(a-1) * (1-err2)) + (err2*(b-1)));
        %result(j)=sign(e(2,j))*abs(d21*d22)*(gamma(a(i)+b(i))/(gamma(a(i))*gamma(b(i))));
        result(j)=sign(e(2,j))*abs(d21*d22)*(gamma(a+b)/(gamma(a)*gamma(b)));
    end
    hold on
    plot(e(2,:),result,colors{i,1});
    [v,I]=max(result);
    hold on
    text(e(2,I)+0.01,result(I)+0.5,strcat(num2str(a(i)),',',num2str(b(i))));
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Beta Hebbian Learning%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
innerloop=10;
lratestep=lRate/iters/innerloop;
w = randn(m,col).*0.001; %1xdim

for times=1:innerloop
    for k=1:iters
        x=data(ceil(rand(1,1)*size(data,1)),:)'; %dimx1
        y = w*x;    %1xdim * dimx1 = 1x1
        if y>50
            aa=1
        end
        e=x-(w'*y);     %dimx1 - (dimx1 * 1x1) = dimx1       
        
        %%%Learning rule update%%%
        er=abs(e);
        
        dr1=((er).^(a-2)).*((1-er).^(b-2));dr2=((-(a-1) * (1-er)) + (er*(b-1)));%avoid imaginary abs(1-er)
        update=lRate * y * (sign(e).*abs(dr1.*dr2))'/(gamma(a)*gamma(b)/gamma(a+b));
        %update(er<0.0001)=0;
        w=w+update;
        lRate=lRate-lratestep;
    end  
end
newdata=(data * w');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                SCATTER PLOT                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Scatter plot'); %Figura donde se ploteará el scatter plot

%Genera el scatter plot enfrentando las m primeras componentes (columnas)
%del conjunto de datos entre si
for i=1:m
    for j=1:m 
       if i==j %La diagonal principal del Scatter plot no se plotea
           hold on;
           subplot(m,m,(m*(i-1))+j)
       else
           for k=1:size(classes,1)%plotea cada clase con un color
               hold on;
               subplot(m,m,(m*(i-1))+j), plot(newdata(tnum==classes(k),i),newdata(tnum==classes(k),j),colors{k,1});
               min1=min(min(newdata(:,i))); %Rescala los ejes al minimo y maximo en el eje X e Y
               maxi1=max(max(newdata(:,i)));
               min2=min(min(newdata(:,j)));
               maxi2=max(max(newdata(:,j)));
               axis([min1,maxi1,min2,maxi2]);
           end
       end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           UNA PROYECCION                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotea dos componentes XComp para el eje X, e YComp para el eje Y

x=input('Seleccione la componente del eje x:  ');
y=input('Seleccione la componente del eje y:  ');
XComp=x;
YComp=y;
    
    figure();

    for k=1:size(classes,1)%cada clase con un color distinto
        hold on;
        plot(newdata(tnum==classes(k),XComp),newdata(tnum==classes(k),YComp),colors{k,1});
    end
    title('2D BHL colors')
    
    figure();
    plot(newdata(:,XComp),newdata(:,YComp),'w.');
    for k=1:size(newdata,1)
        text((newdata(k,XComp)+0.001),newdata(k,YComp),num2str(id(k)));
    end
    title('2D BHL')
    
    
%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                3D                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotea en 3D usando las columnas XComp, YComp, ZComp
XComp=1;
YComp=2;
ZComp=3;

figure();
    for k=1:size(classes,1)%cada clase con un color distinto
        hold on;
        plot3(newdata(tnum==classes(k),XComp),newdata(tnum==classes(k),YComp),newdata(tnum==classes(k),ZComp),colors{k,1});
    end
    title('3D BHL colors')

figure();
    plot3(newdata(:,XComp),newdata(:,YComp),newdata(:,ZComp),'w.');
    for k=1:size(newdata,1)
        text((newdata(k,XComp)+0.001),newdata(k,YComp),newdata(k,ZComp),num2str(id(k)));
    end
    title('3D BHL')

% scatter3(newdata(:,XComp),newdata(:,YComp),newdata(:,ZComp),S,'filled');
% title('3D BHL')