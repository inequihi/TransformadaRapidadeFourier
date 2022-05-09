% MUESTREAR SEÑAL 
[n , f_n] = MuestrearSenial();
size(n,2);

% HACER FFT CON MATLAB
Matlab_F = Matlab_FFT(f_n, size(n,2));

% LEER txt CON FFT de C++
Mi_FFT();

% GRAFICAR FFTs


function [nT , xn] = MuestrearSenial()
    F = 5; %Frecuencia de entrada
    Fs = 40; %Frecuencia de muestreo
    f = F/Fs; 
    A = 2;  %Amplitud
    Fase = 0;  %theta
    T = 1/Fs;  %Tiempo de muestreo

    %La frecuencia relativa debe estar entre -1/2 y 1/2 asi que veifico eso
    while f>(1/2)
        f = f - 1; 
    end

    %Creacion de la señal continua
    t=(0:.001:1); %Arreglo del tiempo continuo. 1001 muestras
    xt = A*cos(2*pi*F*t+Fase); %Funcion continua

    %Creacion de la señal discreta
    n = (0:Fs); % Cantidad de muestras 
    xn = A*cos(2*pi*f*n+Fase);
    nT = n*T;
    
    %Superpongo graficas
    figure(1)
    plot(t,xt,'LineWidth',1);
    hold on
    stem(n*T,xn,'r');
    grid on; 

end
 
function Matlab_F = Matlab_FFT(f_n, cantMuestras)
     Matlab_F = fft(f_n,cantMuestras)
     figure(2)
     plot(f_n,abs(Matlab_F));
end

function Mi_FFT()
   T = readtable('Xfourier.txt','ReadVariableNames',false);
   
   
end
