function h = plot_3di(time, freq, amp, scale, fig_title)
% plot_3d(time, freq, amp, scale, title)
%
% plot 3d-like picture of a signal
% scale can be 'auto' or real number
%
% slazav, feb 2012

    %% plot 3d-like frequency picture 

    amp=abs(amp);
    amin=min(min(amp));

    if nargin > 3
        if strcmp(scale, 'log')      amp=log(amp);
        elseif strcmp(scale, 'sqrt') amp=sqrt(amp-amin);
        else                         amp=amp-amin;
        end
    end
    hold on
    h = image(freq, time, rb_hot(amp));
    if nargin > 4
        ht=title(fig_title);
        set(ht,'Interpreter','none');
    end
end

function col = rb(a)
    m1 = min(min(a));
    m2 = max(max(a));
    for i=1:length(a(:,1))
      for j=1:length(a(i,:))
        v=6*(a(i,j)-m1)/(m2-m1);
        vn=floor(v);
        dv=mod(v,1);
        switch vn
          case 0; col(j,i,:) = [1-dv,1-dv,   1];
          case 1; col(j,i,:) = [   0,  dv,   1];
          case 2; col(j,i,:) = [   0,   1,1-dv];
          case 3; col(j,i,:) = [  dv,   1,   0];
          case 4; col(j,i,:) = [   1,1-dv,   0];
          case 5; col(j,i,:) = [   1,   0,  dv];
          otherwise
            if v<0; col(j,i,:) = [1,1,1];
            else    col(j,i,:) = [1,0,1];
            end
        end
      end
    end
end

function col = rb_hot(a)
    m1 = min(min(a));
    m2 = max(max(a));
    v=3*(a'-m1)/(m2-m1);
    n = size(a,2);
    col = zeros(size(a,2),size(a,1),3);
    col(:,:,1) = min(v,1); % R
    col(:,:,2) = min(max(v-1,0),1); % G
    col(:,:,3) = min(max(v-2,0),1); % B
end


function col = rb_g(a)
    s = 1/max(max(amp));
    col(:,:,1) = 1-amp';
    col(:,:,2) = 1-amp';
    col(:,:,3) = 1-amp';
end


