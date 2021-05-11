function [ Gz ] = BatchLS (input,output, lb,la ,d)

na=la-1;
nb=lb-1;

%% outputs:
% Gz              : discreate transfer function

%% Function body
for i=1:length(input)
    for j=1:na
        if i-j <=0
            phiT(i,j)=0;
        else
            phiT(i,j)=[-output(i-j)];
        end
    end
    for j=0:nb
        if i-j-d <= 0
            phiT(i,j+1+na)=0;
        else
            phiT(i,j+1+na)=[input(i-j-d)];
        end
    end
end

Ts=0.2;

Theta_hat=inv(phiT'*phiT)*phiT'*output';
Gz=tf([Theta_hat(na+1:end)'],[1,Theta_hat(1:na)'],Ts);

    y1=simulator(input,[zeros(1,d),Theta_hat(na+1:end)'],[1,Theta_hat(1:na)']);
    
    figure;
    subplot(3,1,1:2)
    set(gcf,'color','w')
    plot(0:length(output)-1,output,0:length(y1)-1,y1,'-o','linewidth',2);
    grid on;
    ylabel('y, y_e_s_t','fontsize',18);
    legend('y','y_e_s_t')
    subplot(3,1,3)
    plot(0:length(y1)-1,abs(output-y1))
    grid on;
    xlabel('samples','fontsize',18);
    ylabel('y-y_e_s_t','fontsize',18);
    legend('error')

end