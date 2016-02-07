% %Calculation for heterogeneous ice nucleation rate
%
%close all
function KnopfHeadFunction(inputFile, A)
    clc

    format short

    %A=load('100RH_LEO.dat');
    obj=A(:,1); % number of freezing event
    aw=A(:,2); %water activity
    daw=A(:,3); %water activity shift, Koop et al. 2000
    Tf=A(:,4); % freezing temperature
    Tm=A(:,5);  %melting temperature
    fr=A(:,6);  %cooling rate, to get time
    %mr=A(:,9); %melting rate
    sa=A(:,7);  % surface area of IN in droplet
    sa_mean = mean(sa);

    Tmelt=mean(Tm);
    nuc=length(obj);
    Tlast=min(Tf);
    p=0;
    fr=fr./60;
    fr_mean=mean(fr);
    aw=mean(aw);
    sa_mean=mean(sa);

    %%
    dT=-.2; %temperature interval, depends on T uncertainty
    x2=0;x1=0;
    %sumx=0;
    for T=Tmelt:dT:Tlast
        p=p+1;
        x1=nuc;x2=nuc;%number of nucleation event
        sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
        for i=1:nuc
            x1=x1-(T<Tf(i)); %number of liquid samples at the start of the T interval
            x2=x2-((T+dT)<Tf(i)); %number of liquid samples at the end of the T interval
            sum1=sum1+(T>Tf(i))*((T+dT)<=Tf(i))*(T-Tf(i))/fr(i); %inside the interval
            sum2=sum2+((T+dT)>Tf(i))*sa(i); %end of interval
            sum3=sum3+(T>Tf(i))*((T+dT)<=Tf(i))*(T-Tf(i))*sa(i)/fr(i); %inside the interval
            sum4=sum4+(T>Tf(i))*((T+dT)<=Tf(i))*Tf(i); %inside the interval
            sum5=sum5+((T+dT)>Tf(i))*-dT/fr(i); %end of interval
        end
        x=x1-x2; %number of frozen samples at end of interval
        xlist(p,:)=x;
        %sumx=sumx+x; %check to see if the frozen events add up.
        t_tot=sum5*x2+sum1;
        a_t=-dT/fr_mean*sum2+sum3;
        % Lior Addition start
        a = a_t/t_tot;
        n(p,:) = x/a; % excluding the time addition
        % Lior Addition end
        Temp(p,:)=sum4/x;
        w(p,:)=x/t_tot;
        j(p,:)=x/a_t;
        jcalc(p,:)=x/a_t;
    end
    ffcalc=1-exp(-mean(sa)/mean(fr).*cumsum(-jcalc.*dT));

    %This deletes integration values which are inf, not real, or =0.
    lnw=log(w);det=isfinite(lnw);q=0;
    for i=1:length(Temp)
        if det(i)==0
            Temp(i-q)=[];
            xlist(i-q)=[];
            lnw(i-q)=[];
            w(i-q)=[];
            j(i-q)=[];
            n(i-q)=[];
            ffcalc(i-q)=[];
            q=q+1;
        end
    end
  
    Temp=round(Temp.*10)./10;
    
    %% calculate Ns - Lior Segev addition based on D Osullivan 2014, Murray 2011, Niemand 2012
    % ns = -1/sa_mean * ln ( 1 - ni(T)/ntot )
    ns = -1/sa_mean * log(1 - ffcalc);
    
    H = [Temp, j, ns, ffcalc];
    
    figure(2)
    plot(Temp, ffcalc);
    
    %Figures
    figure(1)
    subplot(1,3,1)
    semilogy(Temp,j,'o')
    xlabel('Temp (K^{\circ})');
    ylabel('J_{het}  (cm^{-2}\cdot s^{-1})');
    hold off
    title('J_{het}')
    subplot(1,3,2)
    semilogy(Temp,w,'o')
    hold off
    title('\omega_{het}')
    xlabel('Temp (K^{\circ})');
    ylabel('\omega_{het} (s^{-1})');
    subplot(1,3,3)
    semilogy(Temp,ns,'o')
    hold off
    title('n_{s}')
    xlabel('Temp (K^{\circ})');
    ylabel('n_{s} (cm^{-2})');
    
    %% write to results to file
    % Make file name
    [pathstr,name,ext] = fileparts(inputFile);
    filenameOutput = strcat(name,'-output.dat');
    fileID = fopen(filenameOutput,'w');
    fprintf(fileID,'%s\n', 'Temp(Kelvin), Jhet(1/(cm^2*s)), ns(1/cm^2), FrozenFraction');
    fclose(fileID);
    dlmwrite(filenameOutput,H,'-append','precision','%12.5e')

end

