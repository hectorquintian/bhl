function plotdata(newdata,data,var,XComp,YComp)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           UNA PROYECCION                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotea dos componentes XComp para el eje X, e YComp para el eje Y

figure();
        plot(newdata(:,XComp),newdata(:,YComp),'.b');
    title('2D colors')
    
    figure();
    plot(newdata(:,XComp),newdata(:,YComp),'w.');
    for k=1:size(newdata,1)
        text((newdata(k,XComp)+0.001),newdata(k,YComp),num2str(data(k,var)));
    end
    title('2D')

end

