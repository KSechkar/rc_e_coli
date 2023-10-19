
load exp_selfact_oscill.mat

for i=1:length(LB)
    i
    subplot(3,1,1)
    hold on
    box on
    errorbar(time,OD_Lara(:,i),OD_Lara_err(:,i))

    subplot(3,1,2)
    hold on
    box on
    errorbar(time,GFP_Lara(:,i)./GFP_Lara(1,i),GFP_Lara_err(:,i)./GFP_Lara(1,i))

    subplot(3,1,3)
    hold on
    box on
    plot(time(1:end-1),(OD_Lara(2:end,i)-OD_Lara(1:end-1,i))./(OD_Lara(1:end-1,i).*(time(2:end)-time(1:end-1))))
end

subplot(3,1,1) 
xlabel('time (h)')
ylabel('OD') 

legend([repmat('LB%=',9,1) num2str(LB') repmat('%',9,1) ])

subplot(3,1,2) 
xlabel('time (h)')
ylabel('GFP/OD')