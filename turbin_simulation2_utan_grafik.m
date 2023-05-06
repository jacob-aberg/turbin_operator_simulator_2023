
  
 
function [T,t] = Function(v,F)


    plotta = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parametrar
    dt = 0.1; % s tidsteg
    simulation_tid = 120; % s simulationsstid

% Broms - parameterar
    F_broms = 0; % N bromskraft
    bromsskiva_radie = 2*0.5; % m (effektiv radie skivbroms)
    bromsskiva_densitet = 7850; % kg / m^3
    bromsskiva_area = 4*bromsskiva_radie^2*pi; % m^2
    bromsskiva_tjocklek = 0.05; % m
    bromsskiva_massa = bromsskiva_area * bromsskiva_tjocklek * bromsskiva_densitet; % kg 
    
    spec_varme_kap = 420; % J/kg*K (specific värme kapacitet för bromsskivans material)
    konvektion_koefficent = 100; % W/m^2*K (värmeledningskoefficient / convective heat coefficient)
    start_temp = 20; % C (start temp)
    omgivning_temp = 20; % C (omgivnings temp)
    gear_box_ratio = 100;
    
    my = @(T)   0.55                        .* (T <= 200) + ...
                ( 0.675 - 6.25*10^-4 .* T ) .* (T >  200 & T <= 600) +...
                0.3                         .* (T >  600 & T <= 1200)+...
                0                           .* (T >  1200) ;


% Vindkraftverkets specs
    blad_lngd = 50; % Blade length (m)
    blad_massa = 0.22 * blad_lngd *10^3;
    I_blad = 1/3 * blad_massa * blad_lngd^2 ;%*0.5;
    I = 3*I_blad;
    svept_Area = pi * blad_lngd^2;
    luftdensitet = 1.225; % kg / m^3
    turbin_friktion = 1.1*10^7;
    turbolens = 0.2; % variationen i vindhastighet per tidsteg
    
    % skapar uppskattade effekt-coefficent-modeller
        Cp0 = [0.1, 0.15, 0.25, 0.37, 0.42 , 0.45, 0.4 0.3, 0.21 , 0.15 ];
        CP = interp1( 1:1:numel(Cp0) , Cp0, linspace( 1, numel(Cp0), numel(Cp0)*10 ), 'linear'); %       ..-*¨¨¨-.
    
    % funktion för Cp
        Cp_kurva = @(lambda) diskret_funktion(CP,lambda,10,0);
    % funktion för vindkraftseffekt-ekvationen ( Effekt = RörelseEnergi * flöde * Cp )  
        turbin_effekt = @(v,Cp) 0.5 .* luftdensitet .* Cp .* svept_Area .* v.^3;


% Intitierar vektorer med fysikaliska data
    tid = 0:dt:simulation_tid;  
    broms_moment = zeros(1, length(tid));
    temperatur = zeros(1, length(tid));
    heat_map = zeros(100, 1);
    Ang_momentum = zeros(size(tid));
    w = zeros(size(tid));
    vind_hastighet = zeros(size(tid));
    vind_Effekt = zeros(size(tid));
    vind_moment = zeros(size(tid));
    total_moment = zeros(size(tid));
    friktionstal = zeros(size(tid));
    TSR_vec = zeros(size(tid));
    Cp_vec = zeros(size(tid));


% begynnelsevärden


 F = F ;%10^6 * input('Ange bromskraft [MN]:>> ');
 vind_hastighet(1) = v; %input('Ange vindstyrka [m/s]: >> ');


    Tip_speed_ratio = 7;
    w(1) = 4.670;
    Cp = Cp_kurva(Tip_speed_ratio);
    Ang_momentum(1) = 1.2843*10^8;
    temperatur(1) = start_temp;
    friktionstal(1) = my(temperatur(1)); 
    theta = 0;
    vind_Effekt(1) = turbin_effekt(vind_hastighet(1),Cp);
    vind_moment(1) = blad_lngd * vind_Effekt(1) / vind_hastighet(1);

    broms_klocka = 0;
    turbolens = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
for i = 2:length(tid)

% Vind
    vind_hastighet(i) = vind_hastighet(i-1) + turbolens * sqrt(dt) * randn ;

    Tip_speed_ratio = 7;%TSR( vind_hastighet(i-1), w(i-1), blad_lngd)
    Cp = Cp_kurva(Tip_speed_ratio);
    TSR_vec(i) = Tip_speed_ratio;
    Cp_vec(i) = Cp;

    vind_Effekt(i) = 0.5 * luftdensitet * Cp * svept_Area * vind_hastighet(i)^3;
    vind_moment(i) = 0.5 * luftdensitet * Cp * svept_Area * vind_hastighet(i)^2 * blad_lngd;
    %w_guess = Tip_speed_ratio * vind_hastighet(i) / blad_lngd;
    %if w_guess ~=0 
    %vind_moment(i) = vind_Effekt(i) / w_guess ;
    %end

% bromsen

    if tid(i) > 3
        F_broms = F;
    end

    % beräknar μ för nuvarande temperatur
    friktionstal(i) = my(temperatur(i-1)) ;
        if w(i-1)==0
            % statisk friktion
            friktionstal(i) = friktionstal(i) * 1.45;
        end


    % Beräknar bromsmomentet
    broms_moment(i) = friktionstal(i) * F_broms * bromsskiva_radie  *-1*sign(vind_moment(i));
    % Begränsar bromsmomentet så det inte kan överstiga vindmomentet ,+ motverkande riktining (bromsmoment är reaktivt)
    if Ang_momentum(i-1) <= 0
    broms_moment(i) = min( abs(broms_moment(i)) , abs( vind_moment(i) ) ) *-1*sign(vind_moment(i));
    end
% Bromsens temperatur

    % Värmeutveckling enligt E = ΔΤ * m * C
    heat_generated = abs(broms_moment(i))  * gear_box_ratio*(abs(w(i-1)) * dt) / ( spec_varme_kap * bromsskiva_massa);

    % Värmeledning via konvektion E = h * A * (T1-T2)
    E_konvektion = konvektion_koefficent * bromsskiva_area * ( temperatur(i-1) - omgivning_temp ) ;
    heat_lost = E_konvektion * dt / (spec_varme_kap * bromsskiva_massa); 
    
    % Uppdaterar temperatur
    temperatur(i) = temperatur(i-1) + heat_generated - heat_lost;
    
    % Uppdaterar heat map
    heat_map = linspace(temperatur(i), 0, 100)';


% Rörelsemängden   

    % beräknar generatorns/turbins/luftens friktionsmomentet utifrån linjär modell
    friktion_moment = turbin_friktion * abs(w(i-1))  *-1*sign(vind_moment(i));

    % beräknar momentsumman
    total_moment(i) = vind_moment(i) + gear_box_ratio * broms_moment(i) + friktion_moment; 
    
 
    % integrerar momentet i tiden för att få rörelsemängdsmomentet
    Ang_momentum(i) = Ang_momentum(i-1) + total_moment(i)*dt ;

    if Ang_momentum(i) <=0
        Ang_momentum(i) = 0;
    end
    % Beräknar vinkelhastigheten utifrån röreslemängdens och tröghetens moment
    w(i) = Ang_momentum(i)/I;
    

    if w(i) < 10^(-4) || ~isfinite(w(i))
        w(i) = 0;
    end

    if i>4
    if all(w(i-3:i)<0.01 ), Ang_momentum(i)=0; end, end

% klocka
    if F_broms ~= 0 && w(i)~= 0
        broms_klocka = broms_klocka + dt;
    end

end




if Ang_momentum(i) == 0
%disp("broms tid:  "+ string( broms_klocka ) + " s" )
t = broms_klocka;

else
    %disp("lyckades ej bromsa")
    t = Inf;
end

%disp("max temp:  " + string( max(temperatur) )+ " °C" )
T = max(temperatur);




if plotta

if any(Ang_momentum == 0)
    i = round( find(Ang_momentum == 0,1) * 1.2 ) ;
else 
    i = simulation_tid/dt;
end

% Grafer
    main_figure = figure;
    main_figure.Position = [600,200,0.4*[1920,1080]];
    % 
    %heat map
    ax1 = axes('Position', [0.1   0.80    0.8    0.12 ] );
    % bromsmoment
    ax3 = axes('Position', [0.1    0.60    0.35    0.12] );
    ax32 = axes('Position',[0.55   0.60    0.35    0.12 ]);
    % temp
    ax2 = axes('Position', [0.1   0.35    0.35   0.12] );
    % friktion 
    ax4 = axes('Position', [0.552   0.35   0.35   0.12] );
    % vind
    ax5 = axes('Position', [0.1   0.10   0.80   0.15] );

% updaterar alla plots 
    imagesc(ax1,heat_map);
    colormap(ax1,'jet');
    colorbar;
    caxis(ax1,[0, 100]);
    title(ax1,'Bromskiva',FontSize=10);
   
    plot(ax2,tid(1:i), temperatur(1:i),LineWidth=2,Color='r');
    xlabel(ax2,'Tid (s)');
    ylabel(ax2,'Temperatur (°C)');
    s = 'Bromsskivans värmeutveckling';%+ string(round(temperatur(i))) +' (°C)';
    title(ax2,s,FontSize=10);

    plot(ax3,tid(1:i), -1*broms_moment(1:i),LineWidth=2);
    xlabel(ax3,'Tid (s)');
    ylabel(ax3,' Moment (Nm)');
    title(ax3,'Bromsande moment',FontSize=10);

    plot(ax32, tid(1:i), Ang_momentum(1:i) ,LineWidth=2,Color='#4b248c')
    title(ax32,'Rörelsemängdens moment',FontSize=10)
    xlabel(ax32,'Tid (s)');
    s = "Rörelsemängds-"+newline + "moment" + newline + "(kg m^2/s)";
    ylabel(ax32,s);
    
    plot(ax4,tid(1:i),friktionstal(1:i),LineWidth=2);
    title(ax4,'Friktionstal bromsskiva-bromsbelägg ',FontSize=10);%+string(round(friktionstal(i),2)),FontSize=12);
    ylim(ax4,[0,1])
    xlabel(ax4,'Tid (s)')
    ylabel(ax4, 'μ   ',Rotation=0)

    % plot(ax11,tid(1:i),TSR_vec(1:i),LineWidth=2);
    % title(ax11,'Tip-Speed-ratio (λ)',FontSize=10);
    % ylim(ax11,[-0.1,11])
    % plot(ax12,tid(1:i),Cp_vec(1:i),LineWidth=2,Color='#4b248c');
    % title(ax12,'Power coefficent (C_p) ',FontSize=10);
    % ylim(ax12,[0,0.6])


    plot(ax5,tid(1:i),vind_hastighet(1:i),LineWidth=1,DisplayName='Vindhastighet (m/s)' )
    hold(ax5,"on")
    plot(ax5,tid(1:i),10^-6*vind_Effekt(1:i),LineWidth=1,DisplayName='Effekt (MW)',Color='green' )
    plot(ax5, tid(1:i), w(1:i) ,LineWidth=1,DisplayName='ω (rad/s)')

    hold(ax5,'off')
    title(ax5,'Vindhastighet, effekt, och vinkelhastighet',FontSize=10);
    legend(ax5,'Location', 'northwest');
    drawnow
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = diskret_funktion(Y,x,k,lower_limit)
% k är antalet elemenet mellan funktionsvärden för
% heltalsargument i Y-vektorn
% lower_limit är funktionsvärdet då x<1
if x < 1
    y = lower_limit;
elseif x > numel(Y)/k
    y = Y(end);
else
    idx = k * round(x, k/10);
    if idx ~= 0 && isfinite( idx )
        y = Y( idx );
    else
    y = lower_limit;
    end
end
end


function lambda = TSR(v,w,r)
    if v ~= 0
        lambda = w * r / v ;
    else
        lambda = 0; 
    end
end


end