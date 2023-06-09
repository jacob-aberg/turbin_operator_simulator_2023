
clc,clear all, close all


while true
    fprintf('\nVälkommen till vindkraftslobbyn. \n Vad vill du göra?\n')
    fprintf('1. Testa broms med bromskraft och vindstyrka\n')
    fprintf('2. 3D visualisering\n')
    fprintf('3. Avsluta\n')
    choice = input('>> ');
    
    if choice == 1
            while true
        F = 10^6* input('Ange bromskraft [MN] >> ');
        v = input('Ange vindstyrka [m/s] >> ');

        [temp, tid] = turbin_simulation2_utan_grafik(v , F);

        if tid ~=Inf
        disp("broms tid:  "+ string( tid ) + " s" )
        else
        disp("lyckades ej bromsa")
        end
        disp("max temp:  " + string(temp)+ " °C" )
        
        fprintf(['\nigen?\n' ...
            '\n0. nej' ...
            '\n1. ja\n'])
        if ~input('>> ')
        break
        end
            end
    elseif choice == 2
        visualisera;
    
    elseif choice == 3
        % Exit program
        disp('Hejdå')
        break
    else
        disp('Invalid. försök igen')
    end
end

function visualisera()
disp('Ange önskad upplösning, heltal mellan t.ex. 10 och 1000) (påverkar beräkningstid)')
s = input('>> ');

vind_vektor = linspace(1,40,s); % 1:1:40;
kraft_vektor = 10^6*linspace(1,5,s) ;%[ 1:0.1:5 ] ;

treD_matris = zeros(numel(vind_vektor),numel(kraft_vektor),2);

for i = 1:1:length( vind_vektor )
    for j = 1:1:length( kraft_vektor )
        [temp, tid] = turbin_simulation2_utan_grafik(vind_vektor(i) , kraft_vektor(j) );
        treD_matris(i,j,1:2) = [temp , tid];

    end
end

Z = treD_matris;
Z1 = treD_matris(:,:,1)';
Z2 = treD_matris(:,:,2)';

x = linspace(  min(vind_vektor) , max(vind_vektor), size(Z, 1));
y = linspace(  min(kraft_vektor), max(kraft_vektor), size(Z, 2));

[X,Y] = meshgrid(x, y);

main_figure = figure;


subplot(1,2,2)
surf(X,Y,log10(Z2));
clim([0 8])
xlabel('Vindstyrka [m/s]')
ylabel('Bromskraft [N] ')
zlabel('Tid (s) ')
title('Inbromsningstid')

Z2_log = log10(Z2);


clim([min(Z2_log(:)), max(Z2_log(:))])
cb = colorbar;

ticks = log10( [0:1:10] ); 
ticks = log10( [0 0.25 0.5 1 2 4 6 8 10] );
cb.Ticks = ticks;%[0.5 1 2 5 10 20]; % set desired tick locations
cb.TickLabels = string(10.^(ticks)); %[0.5 1 2 5 10 20]); % set corresponding tick labels



colormap(jet); % Set the colormap
ylabel(cb,'Tid (s)')


subplot(1,2,1)
surf(X,Y,Z1);
clim([20 1300])
xlabel('Vindstyrka [m/s]')
ylabel('Bromskraft [N] ')
zlabel('Temperatur °C')
title('Max temperatur bromsskiva')

cb2 = colorbar;
colormap(cb2,jet); % Set the colormap
ylabel(cb2,'Temperatur (°C)' )

end

