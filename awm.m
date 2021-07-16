clear all

waa = 0.076;
chi = 0.014; 
gamma = 0.3;
kappa = 0.058;
tasa = 0; 
lambda = kappa/(1-kappa);
lado = 100;
N = lado^2;
W = 1;
mu = W/N;
shift = lambda*mu;  
agente(1:lado,1:lado)=mu;
bias = 0;

ejex = linspace(1,N,N);
ejex2 = linspace(1/N,1,N);

count=0;

figure(1)

while(true)

  count = count + 1;
  
  %%%%%%%%%%%%%%%%%%%% MAPA %%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,1)    
  imagesc(agente,[-0.00025,30*mu]);colormap(hot);colorbar;
  title('Riqueza de los distintos agentes',"fontsize",16);
  
  %%%%%%%%%%%%%%%%%%%% GINI %%%%%%%%%%%%%%%%%%%%%
  
  lorentz = sort((vec(agente))./W);
  lorentz = cumsum(lorentz);  
  % Find gini
  cruce = find(lorentz>0,1,'first');
  B = trapz(ejex2(cruce:end),lorentz(cruce:end));
  if(cruce>1)
    C = abs(trapz(ejex2(1:(cruce-1)),lorentz(1:(cruce-1))));
  else
    C = 0;
  end
  A = 0.5 - B;
  gini(count) = (A-C)/(A+B-C);
  
  subplot(2,2,3) ## Lorentz Curve
  plot(ejex2,lorentz);
  title('Curva de lorentz',"fontsize",16);
  xlabel('Fraccion acumulada de agentes F(w)',"fontsize",14);
  ylabel('Fraccion acumulada de riqueza L(w)',"fontsize",14);
  grid on;
  axis([0 1 -0.2 1])
##  disp(count);

  subplot(2,2,4)  ## Gini coefficient
  plot(gini)
  title(sprintf('Coeficiente de Gini: %4.3f',gini(count)),"fontsize",16);
  xlabel('Iteraciones de evolucion temporal',"fontsize",14);
  ylabel('Coeficiente de Gini',"fontsize",14);
  
  
  pause(0.001);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%% ALGORITMO %%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for k=1:N %% Hago N transacciones antes de mostrar resultados gráficamente
    
    ind = ceil(lado*rand(1,4));
    i1 = ind(1);
    i2 = ind(2);
    j1 = ind(3);
    j2 = ind(4);
    
    
      if((i1!=i2)||(j1!=j2))

        bias = waa*(agente(i1,j1)-agente(i2,j2))*N/(sqrt(gamma)*W);

        p = (1+bias)/2;
        r = rand;
        
        if(r < p)
          eta=1;
        else
          eta=-1;
        end  
        
        % SHIFT
        agente(i1,j1) = agente(i1,j1) + shift;
        agente(i2,j2) = agente(i2,j2) + shift;
        
        % Cantidad a transar sin impuestos
        cantidad = min(agente(i1,j1),agente(i2,j2));
        
        delta = sqrt(gamma)*cantidad*eta;
        
        % Transaccion
        agente(i1,j1) = agente(i1,j1) + delta;        
        agente(i2,j2) = agente(i2,j2) - delta;
        
        % Impuesto
        agente(i1,j1) = agente(i1,j1) + chi*(mu+shift-agente(i1,j1));        
        agente(i2,j2) = agente(i2,j2) + chi*(mu+shift-agente(i2,j2));
        
        
        % SHIFT BACK
        agente(i1,j1) = agente(i1,j1) - shift;
        agente(i2,j2) = agente(i2,j2) - shift;
        
                
        %% Condicion periodica de contorno
        
        %% Redistribución a vecinos próximos si tiene más de riqueza media
        if(agente(i1,j1)>mu)
          arr = i1-1;
          aba = i1+1;
          izq = j1-1;
          der = j1+1;
          
          if(i1==1)
            arr = lado;
          end
          if(i1==lado)
            aba = 1;
          end
          if(j1==1)
            izq = lado;
          end
          if(j1==lado)
            der = 1;
          end
            
          agente(arr,j1) = agente(arr,j1) + (tasa/4)*agente(i1,j1);
          agente(aba,j1) = agente(aba,j1) + (tasa/4)*agente(i1,j1);
          agente(i1,izq) = agente(i1,izq) + (tasa/4)*agente(i1,j1);
          agente(i1,der) = agente(i1,der) + (tasa/4)*agente(i1,j1);
          agente(i1,j1) = (1-tasa)*agente(i1,j1);
          
        end
        if(agente(i2,j2)>mu)
          arr = i2-1;
          aba = i2+1;
          izq = j2-1;
          der = j2+1;
          
          if(i2==1)
            arr = lado;
          end
          if(i2==lado)
            aba = 1;
          end
          if(j2==1)
            izq = lado;
          end
          if(j2==lado)
            der = 1;
          end
            
          agente(arr,j2) = agente(arr,j2) + (tasa/4)*agente(i2,j2);
          agente(aba,j2) = agente(aba,j2) + (tasa/4)*agente(i2,j2);
          agente(i2,izq) = agente(i2,izq) + (tasa/4)*agente(i2,j2);
          agente(i2,der) = agente(i2,der) + (tasa/4)*agente(i2,j2);
          agente(i2,j2) = (1-tasa)*agente(i2,j2);
        end
      
      end      
  end 
end


## ESTA SECCION NO SE LLEGA A EJECUTAR,  ##
## HAY QUE COPIARLA Y PEGARLA EN EL TER- ##
## MINAL CUANDO PAREMOS EL PROGRAMA PARA ##
## QUE CREE LA GRÁFICA DEL HISTOGRAMA Y  ##
## EL ARCHIVO CON LOS DATOS DE ESTE PARA ##
## SACAR LA LEY DE POTENCIAS             ##


subplot(2,2,2)
auxhist = vec(agente);
auxhist = auxhist./W;
paso = 0.001/20;
for count=1:20
  aux2 = auxhist(auxhist<=(paso*count));
  aux2 = aux2(aux2>(paso*(count-1)));
  y(count) = size(aux2,1);
  a(count) = paso*count;
end
bar(a,y);
title('Histograma de riqueza',"fontsize",16);
xlabel('Fraccion de riqueza',"fontsize",14);
ylabel('Numero de agentes',"fontsize",14);
res(:,1) = a;
res(:,2)=y;
save res.txt res
 