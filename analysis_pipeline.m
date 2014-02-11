function out_str= analysis_pipeline(g_max,r_max,times,voxel_deg)
    %plot_autocorrelation(voxel_deg,g_max); title('Autocorrelation Green');
    %return;
    %%%%
    %%%%
    %%%%
    %%%%
    %plot green heatmap (red and green)
    out_str.times=times;
    out_str.voxel_deg = voxel_deg;
    out_str.g_max = g_max;
    out_str.r_max = r_max;
    %figure;contourf(times,voxel_deg,g_max);
    %xlabel('time'),ylabel('degrees');title('Green');
    %figure;contourf(times,voxel_deg,r_max);
    %xlabel('time'),ylabel('degrees');title('Red');
    %plot FWHM vs. time (green and red)
    [gmu,gsig]=fit_wrapped_gaussian(g_max);
    fprintf('FWHM green: %.3f +/- %.3f degrees\n',mean(2*sqrt(2*log(2))*gsig),std(2*sqrt(2*log(2))*gsig))
    return;
    [rmu,rsig]=fit_wrapped_gaussian(r_max);
    fprintf('FWHM   red: %.3f +/- %.3f degrees\n',mean(2*sqrt(2*log(2))*rsig),std(2*sqrt(2*log(2))*rsig))
    %%%
    out_str.g_fwhm=2*sqrt(2*log(2))*gsig;
    out_str.r_fwhm=2*sqrt(2*log(2))*rsig;
    %figure;plot(times,2*sqrt(2*log(2))*gsig,'-g',times,2*sqrt(2*log(2))*rsig,'-r');
    %title('FWHM');xlabel('times');ylabel('Width (degrees)');
    %legend('Green','Red','Location','best');
    %angular distance between Green and Red peaks, vs. time
    out_str.tracking = min([abs(gmu-rmu) ; abs(rmu-gmu-360) ; abs(gmu-rmu-360)]);
    fprintf('Tracking  : %.3f  +/- %.3f degrees\n',mean(out_str.tracking),std(out_str.tracking));
    %figure;plot(times,min([abs(gmu-rmu) ; abs(rmu-gmu-360) ; abs(gmu-rmu-360)]),'-b')%,...
    %    times,gmu,'--g',times,rmu,'--r');
    %title('Tracking Accuracy');xlabel('times');ylabel('degrees');
    %legend('abs(diff)','green','red','Location','best');
    return;
    %%%
    %autocorrelation heatmap vs. time (red and green)
    [out_str.g_autocorr_x,out_str.g_autocorr_y]=find_autocorrelation(voxel_deg,g_max); 
    fprintf('Mean correlation distance (green): %g deg\n',mean((out_str.g_autocorr_x * out_str.g_autocorr_y)./sum(out_str.g_autocorr_y)));
    %figure;contourf(data.times,data.g_autocorr_x,data.g_acorr_y);
    %xlabel('time'),ylabel('degrees');title('Autocorrelation Green');
    [out_str.r_autocorr_x,out_str.r_autocorr_y]=find_autocorrelation(voxel_deg,r_max); 
    fprintf('Mean correlation distance   (red): %g deg\n',mean((out_str.r_autocorr_x * out_str.r_autocorr_y)./sum(out_str.r_autocorr_y)));
    %title('Autocorrelation Red');
    %%%
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [xout,yout]=find_autocorrelation(x,y)
        %%%
        for d=1:360
            dat(d,:) = plot_autocorrelation__smoothed_point(d,x,y); %#ok<*AGROW>
            %pause(.1);
        end
        %error('stop');
        %%%
        F=fft(dat);
        C=ifft( F.*conj(F) );
        %C(1:floor(end/2),:);
        xout=1:180;
        yout=C(1:180,:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function val = plot_autocorrelation__smoothed_point(d,x,y)
        if(numel(d)~=1),error('plot_autocorrelation__smoothed_point, d not size 1');end
        ndx=find(x>d,1,'first');
        if(isempty(ndx)),ndx=161;end
        %ndx
        xtmp = circshift(x,-(ndx-1));
        %figure;plot(xtmp)
        %keyboard
        xtmp=xtmp-d;
        xtmp( xtmp<0 ) = xtmp( xtmp<0 )+360;
        %figure;plot(xtmp)
        ytmp = circshift(y,-(ndx-1));
        %surf(times,xtmp,ytmp);
        %plot(xtmp,ytmp(:,1))
        %%%
        weights=zeros(size(xtmp));
        % setup gaussian
        % NOTE: FWHM of 4 deg : sig^2=2/log(2)
        for i=1:length(xtmp)
            if( abs(xtmp(i))>6),break,end
            weights(i) = exp( -(xtmp(i))^2/(4/log(2)));
        end
        for i=length(xtmp):-1:1
            if( abs(360-xtmp(i))>6),break,end
            weights(i) = exp( -(360-xtmp(i))^2/(4/log(2)));
        end
        weights=weights./sum(weights);
        %%%
        for i=1:size(y,2) %for each time point
            val(i) = sum( ytmp(:,i).*weights );
        end
        %%%
        %plot(xtmp,ytmp(:,1),'-b')
        %hold on;
        %plot(xtmp, 255*exp( -(xtmp).^2./(4/log(2)) ),'--r');
        %plot(xtmp, 255*exp( -(xtmp-360).^2./(4/log(2)) ),'--r');
        %plot(0,val(1),'Xk','LineWidth',4)
        %hold off
        %drawnow;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [mu,sig]=fit_wrapped_gaussian(data)
        for x=1:size(data,2)
            mui=fit_wrapped_gaussian__find_mu(voxel_deg,data(:,x));
            %size(data)
            %error('stop');
            %N=sum(data(:,x));
            %z_bar = 1/N*sum(  data(:,x).*exp(1i*voxel_deg/360*2*pi) );
            %mu(x)=angle(z_bar)*360/(2*pi);
            %figure(1);plot(voxel_deg,data(:,x),'-b');
            yvals = data(:,x);
            yvals = yvals-min(yvals);
            yvals=yvals/max(yvals);
            xvals = voxel_deg;
%             z_bar=0;
%             N=0;
%             for i=1:size(data,1)
%                 if(yvals(i)<.75),continue;end
%                 z_bar=z_bar+yvals(i)*exp(1i*xvals(i)/360*2*pi);
%                 N=N+yvals(i);
%             end
%             z_bar=z_bar/N;
%             mu(x)=mod(angle(z_bar)*360/(2*pi)+360,360);
%             
%             %%%
%             mui=find(xvals>=mu(x),1,'first');
            mu(x) = voxel_deg(mui);
            mdi = floor(length(xvals)/2);
            yvals = circshift(yvals,mdi-mui);
            xvals = circshift(voxel_deg,mdi-mui);
            xvals = xvals-xvals(1);
            xvals( xvals<0 ) = xvals(xvals<0) + 360;
            xvals( xvals>360 ) = xvals(xvals>360) - 360;
            %%%
            %fprintf('N=%g mu=%g mui=%g mdi=%g voxel_deg(mui)=%g\n',N,mu(x),mui,mdi,xvals(80));
            %figure(2);plot(xvals,yvals,'-b');
            %%%
            d1 = xvals( find(yvals>.5,1,'first'));
            d2 = xvals( find(yvals>.5,1,'last'));
            sig1=(d2-d1)/(2*sqrt(2*log(2)));
            %y2=exp( -(x-mu).^2 / (sig1^2));
            sigma=sig1^2;
            %%%
            %hold on;plot(xvals,exp( (-(xvals-xvals(80)).^2)/(2*sigma)),'-r');hold off
            %%%
            low = mean( yvals( [1:ceil(length(xvals)/4) 3*ceil(length(xvals)/4):length(xvals) ]  ) );
            i1 = trapz(xvals,yvals);
            %%%
            intlow =  trapz(xvals,ones(size(xvals))*low);
            cntr=1;
            while(intlow>i1)
                low = mean(yvals)-cntr*var(yvals);
                intlow =  trapz(xvals,ones(size(xvals))*low);
                cntr=cntr+1;
            end
            %%%
            y2 = (1-low)*exp( -( xvals - xvals(80) ).^2 / (2*sigma))+low;
            i2 = trapz(xvals,y2);
            i20=i2;
            %%%
            lpcnt=0;
            if(isnan(i1)),error('integral of data is NaN');end
            if(isnan(i2)),sigma=180;end
            while(isnan(i2)|| abs(i1-i2)/i2 > .05 )
                lpcnt=lpcnt+1;
                if(lpcnt==1000)
                    figure;hist(yvals);title(sprintf('mean=%g var=%g',mean(yvals),std(yvals)));
                    figure;plot(xvals,yvals,'-b',xvals,(1-low)*exp( -( xvals - xvals(80) ).^2 / (2*sig1^2))+low,'-r');
                    title(sprintf('i1=%g i2=%g sigma=%g sigma0=%g',i1,i20,sigma,sig1^2));
                    fprintf('low=%g int(low)=%g\n',low,trapz(xvals,ones(size(xvals))*low));
                end
                if(lpcnt>1000)
                    figure;plot(xvals,yvals,'-b',xvals,y2,'-r');
                    title(sprintf('i1=%g i2=%g sigma=%g sigma0=%g',i1,i2,sigma,sig1^2));
                    error('infinite loop');
                end
                if(i1>i2)
                    sigma = sigma*2;
                else
                    sigma = sigma*0.9;
                end
                %low = mean( y( [find(x<mu-3*sqrt(sigma)) ;find(x>mu+3*sqrt(sigma))]  ) );
                y2 = (1-low)*exp( -( xvals - xvals(80) ).^2 / (2*sigma))+low;
                i2 = trapz(xvals,y2);
                %figure(3);plot(xvals,yvals,'-b',xvals,y2,'-r');
                %title(sprintf('i1=%g i2=%g',i1,i2));
                %pause(1);
            end
            %figure(3);plot(xvals,yvals,'-b',xvals,y2,'-r');
            %title(sprintf('i1=%g i2=%g',i1,i2));
            %error('stop');
            sig(x)=sqrt(sigma);
        end
    end
    function mui=fit_wrapped_gaussian__find_mu(x,y)
        szy=length(y);
        dat = [fliplr(y') y' y'];
        %figure(10);plot(x,y,'-b');
        %for itr=1:7
            itr=3;
           winsz=2.5*2^itr;
           cyt = conv(dat,gausswin(szy,winsz));
           cytsz = length(cyt);
           start=floor(cytsz/2)-floor(szy/2);
           cy = cyt(start:start+szy-1);
           cy=255*cy/max(cy);
           %plot(x,y,'-b',x,cy,'--r');title(sprintf('itr=%g',itr));
        %   pause(1);
        %end
        mui=find(cy==max(cy),1,'first');
        %hold on;plot(x(mui),cy(mui),'ob');hold off
        return 
    end

    %%%%%%%%%
%     function [mu,sig]=fit_wrapped_gaussian(data)
%         for x=1:size(data,2)
%             N=sum(data(:,x));
%             z_bar = 1/N*sum(  data(:,x).*exp(1i*voxel_deg/360*2*pi) );
% %             z_bar = 1/N*sum(  data(:,x).*exp(1i*voxel_deg(x)/360*2*pi) );
%             % mu=arg( z_bar )
%             mu(x)=angle(z_bar)*360/(2*pi);
%             %rbar2=z_bar*conj(z_bar);
%             %r2e=N/(N-1)*(rbar2-1/N);
%             %sig(x) = log(1/r2e)/(2*pi)*360;
%             % var(x) = 1 - E[ cos(x_i - mu)]
%             sig(x) = (1 - 1/N.*sum(data(:,x).*cos(voxel_deg-mu(x)))) /(2*pi)*360;
%         end
%     end
end
