function [BPM, DOM_PEAK, N_prev]= fractal_algorithm(sig,BPM_00,srate,p,idnb, N_prev, N_LMS, N_RLS,see,BPM_orig, BPM_DOM)


   %%
    %taking individual array
    PPG = .5*(clean_up(sig(1,:),srate)+clean_up(sig(2,:),srate)); %adding both channel
    XData = clean_up(sig(3,:),srate);
    YData = clean_up(sig(4,:),srate);
    ZData = clean_up(sig(5,:),srate);
    L0 = length(BPM_00);
    
%     if see==1
%         EPP = [clean_up(sig(1,:),srate); clean_up(sig(2,:),srate); PPG];
%         save_str = sprintf('Data Dump 2\\Pre_Process_%02d_%03d', idnb, p);
%         save(save_str, 'EPP');
%         
%         factor = 125/16385*60;
%         NFFT = 16384;
%         XFFT = (0:NFFT-1)*factor;
%         figure('Visible','off');
%         io = round(BPM_orig/factor)+1;
%         
%         FF = abs(fft(EPP(1,:),NFFT)).^2/1000;
%         [Val,MI] = max(FF(1:NFFT/2));
%         subplot(3,1,1),plot(XFFT,FF,'k');
%         hold on, plot(XFFT(io),FF(io),'ob');
%         hold on, plot(XFFT(MI),FF(MI),'dr');
%         xlim([0 200]);
%         ylabel('(a)');
%         title('Periodogram of PPG Signal 1 after Filtering');
%         
%         FF = abs(fft(EPP(2,:),NFFT)).^2/1000;
%         [Val,MI] = max(FF(1:NFFT/2));
%         subplot(3,1,2),plot(XFFT,FF,'k');
%         hold on, plot(XFFT(io),FF(io),'ob');
%         hold on, plot(XFFT(MI),FF(MI),'dr');
%         xlim([0 200]);
%         ylabel('(b)');
%         title('Periodogram of PPG Signal 2 after Filtering');
%         
%         FF = abs(fft(EPP(3,:),NFFT)).^2/1000;
%         [Val,MI] = max(FF(1:NFFT/2));
%         subplot(3,1,3),plot(XFFT,FF,'k');
%         hold on, plot(XFFT(io),FF(io),'ob');
%         hold on, plot(XFFT(MI),FF(MI),'dr');
%         xlim([0 200]);
%         ylabel('(c)');
%         xlabel('Heart Rate (BPM)');
%         title('Periodogram after adding the two PPG Signals');
%         print(save_str,'-dpng');
%     end
    %%
    %Some constants
    M = 4096/(srate*60);
    %%
    %LMS FILTER DATA
    lms_order=26; %length
    mux = 0.0001; %step size
    muy = 0.0001;
    muz = 0.0001;
    
    %%
    %rls filter data
    rls_lambda=0.999;
    rls_order = 55;
    P0 = 10*eye(rls_order);
    
    %%
    %Affine Projection with Recursive Update
    ap_order = 26;
    muapx = 0.008;
    muapy = 0.008;
    muapz = 0.008;
    
    %%
    if L0==0
        Y0 = initial_estimate3(PPG,XData,YData,ZData,srate);
        N_prev = Y0*M+1; 
%         [a,Nprev] = max(abs(fft(PPG,4096)));
%         Nprev = Nprev(1);
%         Y0 = Nprev/M;
    else
        Y0=BPM_00(end);
        %Nprev = round(Y0*M+1);
    end
    
    %%

    
    %Perfoming LMS
    
    E1 = zeros(1,1000);
    for i = 1:N_LMS
        E1 = E1 + lms_step(PPG,XData,YData,ZData,mux,muy,muz,lms_order);
    end
    E1 = E1/norm(E1); %Normalizing
    

%     %Performing RLS
    
     E2 = zeros(1,1000);
        for i = 1:N_RLS
           E2 = E2 + rls_step(PPG,XData,YData,ZData,rls_lambda,rls_order,P0);
        end
        E2 = E2/norm(E2); %normalizing
     
     %Performing Affine Projection
       E3 = zeros(1,1000);
%      for i = 1:N_RLS
%         E3 = E3 + ap_step(PPG,XData,YData,ZData,muapx,muapy,muapz,ap_order);
%      end
%      E3 = E3/norm(E3); %normalizing

%      if see==1
%         figure,plot(E1);
%         hold on, plot(E2,'r');
%         hold on, plot(E3,'k');
%         drawnow;
%      end

    lambda =  1/2;
    
    %EX = lambda* (E1(1:end-3)+ E2(4:end) + E3(4:end));
    EX = lambda* (E1+ E2 + E3);
%     NFFT = 16384;
%     maxp = 0;
%     maxpp = 0;
%     for i=0:3
%         EX1 = E1(1:end-i) + E2(i+1:end);
%         FFEX1 = max(abs(fft(EX1,NFFT)));
%         
%         if(FFEX1>maxp)
%             maxp = FFEX1;
%             EX = EX1;
%             maxpp = i;
%         end
%         
%         EX1 = E1(i+1:end) + E2(1:end-i);
%         FFEX1 = max(abs(fft(EX1,NFFT)));
%         
%         if(FFEX1>maxp)
%             maxp = FFEX1;
%             EX = EX1;
%             maxpp = -i;
%         end
%     end
    
    
    %EX = E1+E2;
    EX = clean_up(EX,srate);
    [BPM, DOM_PEAK, N_prev]= bpmFft3(EX,4,Y0, srate, N_prev); %%Applying Trivial BPM Calculation
    %[BPM N_prev]= bpmFft4(EX,26,Y0, srate, N_prev); %%Applying Trivial BPM Calculation
    
    
    %if length(BPM_00) == 22
    %    x=1;
    %end
    
%     if see==1
%         EPP = [PPG;XData;YData;ZData;E1;E2;EX];
%         save_str = sprintf('Data Dump\\Comb_%02d_%03d', idnb, p);
%         save(save_str,'EPP');
%         
%         factor = 125/16385*60;
%         NFFT = 16384;
%         XFFT = (0:NFFT-1)*factor;
%         figure('Visible','off');
%         io = round(BPM_orig/factor)+1;
%         
%         FF = abs(fft(PPG,NFFT)).^2;
%         subplot(4,2,1),plot(XFFT,FF);
%         hold on, plot(XFFT(io),FF(io),'o');
%         xlim([0 200]);
%         title('Pre Processed PPG Signal');
%         
%         FF = abs(fft(XData,NFFT)).^2;
%         subplot(4,2,3),plot(XFFT,FF);
%         hold on, plot(XFFT(io),0,'o');
%         xlim([0 200]);
%         title('X Channel Accelerometer Data');
%         
%         FF = abs(fft(YData,NFFT)).^2;
%         subplot(4,2,5),plot(XFFT,FF);
%         hold on, plot(XFFT(io),0,'o');
%         xlim([0 200]);
%         title('Y Channel Accelerometer Data');
%         
%         FF = abs(fft(ZData,NFFT)).^2;
%         subplot(4,2,7),plot(XFFT,FF);
%         hold on, plot(XFFT(io),0,'o');
%         xlim([0 200]);
%         title('Z Channel Accelerometer Data');
%         
%         FF = abs(fft(PPG,NFFT)).^2;
%         subplot(4,2,2),plot(XFFT,FF);
%         hold on, plot(XFFT(io),FF(io),'o');
%         xlim([0 200]);
%         title('Pre Processed PPG Data');
%         
%         FF = abs(fft(EX,NFFT)).^2;
%         subplot(4,2,4),plot(XFFT,FF);
%         hold on, plot(XFFT(io),FF(io),'o');
%         xlim([0 200]);
%         title('Output of Combination of LMS and RLS Filtering');
%         
%         FF = abs(fft(E1,NFFT)).^2;
%         subplot(4,2,6),plot(XFFT,FF);
%         hold on, plot(XFFT(io),FF(io),'o');
%         xlim([0 200]);
%         title('Output of Cascaded LMS filtering');
%         
%         FF = abs(fft(E2,NFFT)).^2;
%         subplot(4,2,8),plot(XFFT,FF);
%         hold on, plot(XFFT(io),FF(io),'o');
%         xlim([0 200]);
%         title('Output of Cascaded RLS filtering');
%         
%         print(save_str,'-dpng');
%         
%         p
%     end
    
    if L0>3
        
        BPM = 0.9*BPM+0.05*BPM_00(end)+0.05*BPM_00(end-1);
        
        %new Code
%         val = 10;
%         if DOM_PEAK-BPM_DOM(end) > val
%             DOM_PEAK = BPM_DOM(end)+val;
%         elseif DOM_PEAK-BPM_DOM(end) < -val
%             DOM_PEAK = BPM_DOM(end)-val;
%         end
%         BPM_ext = BPM_DOM(end-3:end)';
%         xx = [(0:3)' ones(4,1)];
%         a = (xx'*xx)\(xx'*BPM_ext);
%         BPM_est = [4,1]*a;
%         if abs(BPM-BPM_est)>10
%             BPM = DOM_PEAK;
%         end
        %new Code ends
        
%          if((BPM-BPM_00(end))>=5)
%           BPM=BPM_00(end)+5;
%           N_prev  = BPM*M+1;
%          elseif((BPM-BPM_00(end))<=-3)
%           BPM=BPM_00(end)-3;
%           N_prev  = BPM*M+1;
%          end
    end
end