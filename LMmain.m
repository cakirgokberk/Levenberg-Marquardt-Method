% Optimization Theory
% Levenberg Marquardt
% 07/01/2021


 clear; close; clc;

%% Model
a=2                                                                 ;
b=0.7                                                               ;
x = 0.1:0.2:5                                                       ;
e1= 1e-4                                                            ;       % Tolerance
y = log10(a*x)+3*cos(b*x)                                           ;



%% Initial Conditions
a_init = 0.5                                                        ;       %  Initial Parameters 
b_init = 0.5                                                        ;
f_new = (y)' - log10(a_init*x') + cos(b_init*x')                    ;
max_iter = 25                                                       ;
delta = [0;0]                                                       ;
mu = 5                                                              ;
mu_values = [mu]                                                    ;
i_values = [1]                                                      ;



%% Levenberg Marquardt Method
a_new = a_init;
b_new = b_init;

for i = 2:max_iter
    
    J = jacobian(a_new,b_new,x)                                     ;
    f_next = f_new + J*[delta(1,1);delta(2,1)]                      ;
    next_delta = -(pinv(((J')*J + mu*eye(2))))*((J'))*(f_new)       ;
    a_new = a_new - next_delta(1,1)                                 ;
    b_new = b_new - next_delta(2,1)                                 ;
    delta = next_delta                                              ;
    f_prev = f_new                                                  ;
    f_new = (y)' - (log10(a_new*x')+3*cos(b_new*x'))                ;
    Err=((f_new).^2)                                                ;                                             
    
    
     if(sum(f_prev.^2) > sum(f_new.^2))
         
         mu = 0.5*mu                                                ;
         
     else
         
         a_new = a_new + next_delta(1,1)                            ;
         b_new = b_new + next_delta(2,1)                            ;
         f_new = f_prev                                             ;
         mu = 2*mu                                                  ;
         
     end
     
     mu_values  = [mu_values, mu]                                   ;
     i_values   = [i_values, i]                                     ;
     param=[a_new b_new]                                            ;       % Estimated Parameters
     Err_values=Err(1:i)                                            ;       % Errors
     
     
     
    %% Plotting
    figure(1) 
    subplot(3,1,1)
    plot(x,y,'k-.')                             ; 
    hold on                                     ;
    x = 0.1:0.2:5                               ;
    y_new = log10(a_new*x)+3*cos(b_new*x)       ;
    plot(x,y_new,'r')                           ;
    xlabel('X')                                 ;
    ylabel('Function Value')                    ;
    legend('Model Function','Fitted Function')  ;
    title('Fitting')                            ;
    drawnow
    pause(0.007)
    hold off                                    ;
    
    subplot(3,1,2)
    plot(i_values,mu_values, 'b--o')            ;
    xlabel('Number of Iteration')               ;
    ylabel('Value of Mü')                       ;
    title('Value of Mü vs. Iteration')          ;
    hold off                                    ;
    
    subplot(3,1,3)
    plot(Err_values,'m--o')                     ;
    %hold on
    xlabel('Number of Iteration')               ;
    ylabel('Error')                             ;
    title('Error')                              ;
    hold off
    
     if abs(f_new) <e1
     break;
    end
    
    
end
disp('----------------------------------------------------------');
fprintf('Number of Iteration    : %d\n',i)                        ;
fprintf('Model Parameters       : %f     %f \n',a,b)              ;
fprintf('Predicted Parameters   : %f     %f \n',param(1),param(2));
disp('----------------------------------------------------------');



