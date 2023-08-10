clear;

Y25=[3.002849737 8.066909988 22.35681189; 2.615580058 7.794190957 21.67247704; 2.881939562 8.232420866 23.65656795]; % these are N15 rates at 25% of full strength; there are 3 times points and 3 replicate at each time point
Y50=[9.338685229 16.96355479 45.96433913; 9.360884754 20.45955453 49.65061607; 9.832918859	20.4872165	48.74773586]; % these are N15 rates at 50% of full strength; there are 3 times points and 3 replicate at each time point
Y75=[14.67321399 31.05515526 73.98568008; 14.25415075 29.65449586 73.89993511; 14.64291922 29.50779635	70.51195535]; % these are N15 rates at 75% of full strength; there are 3 times points and 3 replicate at each time point
Y100=[24.56148709 48.27908222 110.9907294; 23.31482412 47.76151813 111.6486858; 24.95313797 43.64514589	110.3327729]; % these are N15 rates at 100% of full strength; there are 3 times points and 3 replicate at each time point

% show data to make sure they are correct

disp('input data');
disp(Y25);
disp(Y50);  
disp(Y75);
disp(Y100);

% calculate averages of triplicates

disp('averaged data');
y25=mean(Y25);
y50=mean(Y50);
y75=mean(Y75);
y100=mean(Y100);

disp(y25);
disp(y50);
disp(y75);
disp(y100);


% calculate variance associated with noise in data

v=sum(var(Y25))+sum(var(Y50))+sum(var(Y75))+sum(var(Y100));
v=v/24;
disp('variance of noise');
disp(v);

s=sqrt(v/3); % this is the standard deviation for the average of triplicates

y25store=y25;
y50store=y50;
y75store=y75;
y100store=y100;

muvalues=[];
gvalues=[];

% fit data 1001 times, the first time with no noise
% other times with noise
for ntimes=1:1001;
    if ntimes==1;
        y25=y25store;
        y50=y50store;
        y75=y75store;
        y100=y100store;
    else;
        y25=y25store+normrnd(0,s,1,5);% if there are three time points, add the voise in vectors with three elements.   
        y50=y50store+normrnd(0,s,1,5);
        y75=y75store+normrnd(0,s,1,5);
        y100=y100store+normrnd(0,s,1,5);
    end;

gmat=[0.1:0.01:1.50]; % possible grazing rates per day, I make a larger range.
mumat=[0.1:0.01:1.50]; % possible growth rates per day



SS=10^6; % set SS to very large number initially

for ge=1:length(gmat);
    g=gmat(ge);
    for me=1:length(mumat);
        mu=mumat(me);
        
        % fit 25% of full strenth
        k=mu-0.25*g;
        f75=(exp(k*7.5/24)-1)/k;
        f135=(exp(k*13.5/24)-1)/k;
        f290=(exp(k*29/24)-1)/k;
        a=(f75*y25(1)+f135*y25(2)+f290*y25(3))/(f75^2+f135^2+f290^2);
        yhat75=a*f75;
        yhat135=a*f135;
        yhat290=a*f290;
        yhat=[yhat75;yhat135;yhat290];
        ss=(y25(1)-yhat75)^2+(y25(2)-yhat135)^2+(y25(3)-yhat290)^2;
        
        % fit 50% of full strenth
        k=mu-0.5*g;
        f75=(exp(k*7.5/24)-1)/k;
        f135=(exp(k*13.5/24)-1)/k;
        f290=(exp(k*29/24)-1)/k;
        a=(f75*y50(1)+f135*y50(2)+f290*y50(3))/(f75^2+f135^2+f290^2);
        yhat75=a*f75;
        yhat135=a*f135;
        yhat290=a*f290;
        yhat=[yhat;yhat75;yhat135;yhat290];
        ss=ss+(y50(1)-yhat75)^2+(y50(2)-yhat135)^2+(y50(3)-yhat290)^2;
        
         % fit 75% of full strenth
        k=mu-0.75*g;
        f75=(exp(k*7.5/24)-1)/k;
        f135=(exp(k*13.5/24)-1)/k;
        f290=(exp(k*29/24)-1)/k;
        a=(f75*y75(1)+f135*y75(2)+f290*y75(3))/(f75^2+f135^2+f290^2);
        yhat75=a*f75;
        yhat135=a*f135;
        yhat290=a*f290;
        yhat=[yhat;yhat75;yhat135;yhat290];
        ss=ss+(y75(1)-yhat75)^2+(y75(2)-yhat135)^2+(y75(3)-yhat290)^2;
        
        % fit 100% of full strenth
        k=mu-g;
        f75=(exp(k*7.5/24)-1)/k;
        f135=(exp(k*13.5/24)-1)/k;
        f290=(exp(k*29/24)-1)/k;
        a=(f75*y100(1)+f135*y100(2)+f290*y100(3))/(f75^2+f135^2+f290^2);
        yhat75=a*f75;
        yhat135=a*f135;
        yhat290=a*f290;
        yhat=[yhat;yhat75;yhat135;yhat290];
        ss=ss+(y100(1)-yhat75)^2+(y100(2)-yhat135)^2+(y100(3)-yhat290)^2;
        
        if ss<SS;
            SS=ss;
            gbest=g;
            mubest=mu;
            %disp([SS gbest mubest]);
        end;
    end;
end;

if ntimes==1;
disp('best g');
disp(gbest);
disp('best mu');
disp(mubest);
disp('SS');
disp(SS);
V=SS/6;

F=V/v;
disp('F statistic with 6 and 24 degrees of freedom');
p=1-fcdf(F,6,24);
disp(F);
disp('type I error rate for F');
disp(p); 


plot(y25',yhat(1:3),'k^');
hold on;
plot(y50',yhat(4:6),'r^');
plot(y75',yhat(7:9),'ks');
plot(y100',yhat(10:12),'rs');
xlabel('Y');
ylabel('Calculated Y');
set(gca,'fontsize',12);
hold off;
end;

if ntimes>1;
    gvalues=[gvalues gbest];
    muvalues=[muvalues mubest];
end;
end;

gvalues=sort(gvalues);
muvalues=sort(muvalues);

disp('95% confidence interval for g values')
disp([gvalues(26) gvalues(974)]);

disp('95% confidence interval for mu values')
disp([muvalues(26) muvalues(975)]);

disp(yhat(1:3));
disp(yhat(4:6));
disp(yhat(7:9));
disp(yhat(10:12));
        
        
        

