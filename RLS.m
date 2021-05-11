function Gz = RLS( input, output, lb, la, d, lambda )
na=la-1;
nb=lb-1;
Ts=.2;
nu = na + nb + 1; 

N = length(output); 
n = max(na+1, nb+d+1);



alpha = 10E6;
for i = 1:(n-1)
    P{1,i} = alpha*eye(nu, nu);
end


theta = cell(1, N);
for i = 1:(n-1)
    theta{1, i} = zeros(nu,1);
end


Name = cell(nu, 1);




 for j = n : N
    
        for i = 1 : na 
            if ((j-i)<=0)
                Phi(j, i) = 0;
            else
                Phi(j, i) = -output((j-i));
            end
        end
    
        
        for i = na+1 : nu  
            if ((j-(i-na)-d)<=0) 
                Phi(j, i) = 0;
            else
                Phi(j, i) = input((j-d-i + (na+1)));
            end
        end
        

        phi = Phi(j, :)'; 
        
      
            P{1, j} = 1/lambda*(P{1, j-1} - P{1, j-1}*phi*inv(lambda+phi'*P{1, j-1}*phi)*phi'*P{1, j-1}); 
            K = P{1, j}*phi; % Gain
            theta{1, j} = theta{1, j-1} + K*(output(j)-phi'*theta{1, j-1});
            A = 1;
            for i = 1: na
                A(i+1) = theta{1, j}(i);
            end
            for i = na+1:nu
                B(i-na) = theta{1, j}(i);
            end

        
 end
   
y_estm=simulator(input,[zeros(1,d),B],A);

Theta = cell2mat(theta);
Theta_final = Theta(:, N); 

 
 z = tf('z');
 A = z^(d+nb);
 B = 0;
for i = 1:na
    a(i) = Theta_final(i);
    A = A + a(i)*z^(d-na+nb+na-i);
end
for i = na+1:nu
    b(i-na) = Theta_final(i);
    B = B + b(i-na)*z^(nb - (i-(na+1)));
end

Gz = B/A;
[num, den] = tfdata(Gz);
num = cell2mat(num);
den = cell2mat(den);
Gz = tf(num, den, Ts);

final_time = (N-1)*Ts;
t = 0:Ts:final_time;

figure
subplot(121)
plot(t/Ts, output, 'LineWidth', 2)
hold on
plot(t/Ts, y_estm,'MarkerSize',3, 'Marker','square','color', 'r', 'LineWidth', 2);
xlabel('Samples')
ylabel('Output');
title('Output from Yaw acceleration')
grid on 


subplot(122)
for i = 1 : nu
    hold all
    plot(t/Ts, Theta(i, :), 'LineWidth', 2);
    if (i>=1 & i<=na)
        Name{i, 1} = ['a', num2str(i)]; 
    else
        Name{i, 1} = ['b', num2str(i-na-1)]; 
    end
end
ylim([-10 20])
grid on
xlabel('Samples')
ylabel('Estimated Parameters')
legend(Name)

end
