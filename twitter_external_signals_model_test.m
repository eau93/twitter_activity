function twitter_external_signals_model_test(N,G,k_avg,k_in,k_out,total_time,tau,A,f,E,u)

% testing for github

% this function creates an N x time matrix of ones and zeros, one for
% tweet and zero for no tweets
% to access this matrix type 'load external_signals_tweets.mat'
% the matrix is named external_siganls_tweets_mat
% to open plots for 15_5 division 
% type 'open twitter_external_signals_model_plots_for_15_5_division.fig'
% type 'open csvlist.dat' to get list of connections
% to save the matrix as a .txt type 'load external_signals_tweets.mat'
% 'conten=who;' 'save('external_signals_tweets.txt',conten{:},'-ascii')'
% 'open external_signals_tweets.txt'


% N=100
% G=5
% k_avg=20, case 1: k_in=15 k_out=5, case 2: k_in=19 k_out=1
% time=24
% tau=5
% A=1/2
% f=.1 (fraction of group to get 'kick')
% E=.1
% u=1 (average 'kicks' per day)

% add 24 phase shift
time=total_time+24;

% create a vector of all the users
users=zeros(1,N);
for n=1:N
    users(n)=n;
end

group_size=N/G;

M=connections(N,G,group_size,users,k_in,k_out,k_avg);
map=connections_table(N,G,group_size,k_avg,M);

influence_mat=zeros(time,N);

tweets_mat=zeros(N,time);

indep_tweets_mat=zeros(N,time);

alpha=0;

m=zeros(1,24);
n=zeros(1,24);
a=A*rand;
for i=1:24
    m(i)=(sin((pi()/12)*i)+1);
    if a<m(i)
        n(i)=1;
    else
        n(i)=0;
    end
end
n_s=sum(m);
alpha=u/n_s;

total_kicks_vec=zeros(1,time);
total_kicks_mat=zeros(G,time);
for t=1:time
    
    r_i=A*(sin((pi()/12)*t)+1);
    k=r_i*alpha;
    r=zeros(1,time);
    for i=1:time
        r(i)=A*(sin((pi()/12)*i)+1);
        r_tot=sum(r);
    end
    f_i=(((1/time)*r_tot)/k_avg);
   
    kick=zeros(1,G);
    for g=1:G 
        b=A*rand;
        if b<=k
            kick=1;
        else
            kick=0;
        end
            total_kicks_mat(g,t)=kick;
    end
    if sum(kick)>=1
        kick_time=1;
    else
        kick_time=0;
    end
    total_kicks_vec(t)=kick_time;
    
    [independent_tweets]=independent_activity(N,G,group_size,r_i,t,k,E,f,total_kicks_mat,A);
    
    if t==1
        influence_tweets=zeros(1,N);
    else
        influence_tweets=influence_prob(N,influence_tot,r_i,f_i,k_avg,A);
    end
    
    % total_tweets takes into account the tweets from independent activity
    % and from influenced activity
    total_tweets=independent_tweets+influence_tweets;
    influence_vec=influence(N,k_avg,total_tweets,map);
    influence_mat(t,:)=influence_vec;
    
    % if user tweets decrease influence from previous time step by 1
    for i=1:N
       if t>1
           if total_tweets(i)>0
               total_tweets(i)=1;
           end
           if total_tweets(i)>=1
               if influence_mat(t-1,i)>0
                   influence_mat(t-1,i)=influence_mat(t-1,i)-1;
               else
                   influence_mat(t-1,i)=influence_mat(t-1,i);
               end
           end
       end
    end
    
    % the total influence is the sum of the influences from the previous
    % tau timesteps

    if t<=5
        a=1;
    else
        a=t-tau;
    end

    if t>1
        influence_tot=sum(influence_mat(a:t-1,:),1);
    else
        influence_tot=influence_vec;
    end
    indep_tweets_mat(:,t)=independent_tweets;
    tweets_mat(:,t)=total_tweets;
    
end
days=total_time/24;
total_kicks_vec=total_kicks_vec(25:time)
total_kicks_mat=total_kicks_mat(:,25:time);
k
total_kick=sum(total_kicks_vec)
total_kicks=sum(total_kicks_mat,2)
kicks_per_day=total_kicks/days

tweets_mat;

external_signals_tweets_mat=tweets_mat;

save('external_signals_tweets.mat','external_signals_tweets_mat');
%type('tweets.mat');

indep_tweets_mat;
save('indep_tweets.mat','indep_tweets_mat');
indep_vec=sum(indep_tweets_mat,1);

tweets_mat=tweets_mat(:,25:time);
indep_tweets_mat=indep_tweets_mat(:,25:time);
indep_vec=indep_vec(25:time);

h=gobjects(0);

%make_individual_plots(tweets_mat,indep_vec,time,G,group_size)

figure
make_plots(tweets_mat,indep_vec,total_time,G,group_size,total_kicks_mat,indep_tweets_mat)
h(1)=gcf;
figure
make_plots_minus_seasonality(tweets_mat,indep_vec,total_time,G,group_size,total_kicks_mat)
h(2)=gcf;
figure
make_plots_random_groups(tweets_mat,indep_vec,total_time,G,group_size,N,total_kicks_mat)
h(3)=gcf;

string=sprintf('twitter_individual_model_plots_for_%d_%d_division',k_in,k_out)
saveas(h,string,'fig');
%close all

end

function [independent_tweets]=independent_activity(N,G,group_size,r_i,t,k,E,f,total_kicks_mat,A)
% create a network with the different groups
% these are the conditions at the first time step t=1
independent_tweets=zeros(1,N);
i=0;
y_i=[1:group_size];
fraction=group_size*f;
for g=1:G
    % create a vector of users, inputs will be 1 for tweet or 0 for no
    % tweet
    y=randsample(y_i,fraction);
    kicked_members=sort(y);
    
    kick=total_kicks_mat(g,t);
    
    q=1;
    for j=1:group_size
        %if q>fraction
         %   m=0;
        if j==kicked_members(q)
            m=1;
            q=q+1;
        else
            m=0;
        end
        
        delta=kick*E*m;
        a=rand;
           if g==2
               delta=0;
           end
        if a<r_i+delta
            tweet=1;
        else
            tweet=0;
        end
        
        i=i+1;
        independent_tweets(i)=tweet;
    end
end
independent_tweets;

end

function [M]=connections(N,G,group_size,users,k_in,k_out,k_avg)
% this function connects the users to their inputs

keyType=[];
valueType=[];

map=containers.Map('keyType','int32','valueType','any');

group=zeros(1,group_size);
for i=1:G
    groups{1}(1)=1;
    for k=2:group_size;
        groups{k}(1)=groups{k-1}(1)+group_size;
    end
    for j=2:group_size
        groups{i}(j)=groups{i}(j-1)+1; 
    end
end


for n=1:N
    % determine which group the user is in
    a=ceil(n/group_size);
    
    % find the inputs to the user
    inputs_in=randsample(groups{a},k_in);
    
    % make sure users do not input to themselves
    b=mod(n,group_size);
    if b==0
        b=200;
    end
    
    c=groups{a}(b);
    for i=1:k_in
        if inputs_in(i)==c;
            inputs_in=randsample(groups{a},k_in);
        end
    end
    % find the inputs outside of the group
    users_out=setdiff(users,groups{a});
    inputs_out=randsample(users_out,k_out);
    inputs=[inputs_in inputs_out];

    % map(n) is a vector of all the inputs to user n
    map(n)=inputs;
    
end

vec=zeros(k_avg,1);
use=ones(k_avg,1);
for i=2:N
    vec(1:k_avg)=i;    
    use=[use ;vec];
end
use;

map_mat=map(1);
for i=2:1000
    map_mat=[map_mat map(i)];
end
in=transpose(map_mat);
b=[use in];
csvwrite('csvlist.dat',b);
type csvlist.dat;

filename='csvlist.dat';
M=csvread(filename);

end

function [map]=connections_table(N,G,group_size,k_avg,M)
% this function saves the connections made

keyType=[];
valueType=[];

map=containers.Map('keyType','int32','valueType','any');

inputs=zeros(1,k_avg);
for n=1:N
   for i=1:k_avg
       inputs(i)=M(i+(k_avg*(n-1)),2);
       
   end
   map(n)=inputs;
end

end

function [influence_vec]=influence(N,k_avg,total_tweets,map)
% this function creates a vector with N inputs, the ith component is the number
% of the ith users inputs that tweeted

tweet_vec=zeros(1,k_avg);
influence_vec=zeros(1,N);
for n=1:N
   vec=map(n);
   for i=1:k_avg
      a=vec(i);
      tweet=total_tweets(a);
      tweet_vec(i)=tweet;
   end
   tweet_vec;
   C=sum(tweet_vec);
   influence_vec(n)=C;
end
influence_vec;
end

function [influence_tweets]=influence_prob(N,influence_tot,r_i,f_i,k_avg,A)
% this function determines whether or not a user tweets based on the
% influence from their inputs

influence_tweets=zeros(1,N);
for n=1:N
    a=rand;
    C=influence_tot(n);
    p=1-(1-r_i)*(1-f_i)^C;
    if a<=p
        tweet=1;
    else
        tweet=0;
    end
    influence_tweets(n)=tweet;
end
influence_tweets;
end

function make_individual_plots(tweets_mat,indep_vec,total_time,G,group_size)
% plot of total data over time
figure
% plot the tweets over time in red
vec=sum(tweets_mat,1);
plot(vec,'r')
title('total data: # of tweets over time')
xlabel('time (hours)')
ylabel('number of tweets')
hold on
% plot the independent tweets over time for comparison in blue
plot(indep_vec,'c')
legend('total data','independent tweets','location','NorthEastOutside')

% plot of total data at time t vs total data at time t+1
figure
A_c=vec;
A_c_next=A_c;
A_c_next(1)=[];
A_c(time)=[];
plot(A_c,A_c_next,'o')
title('total data: # of tweets at time t vs # of tweets at time t+1')
xlabel('# tweets at time t')
ylabel('# tweets at time t+1')

% plot the average data
period=time/24;
A=zeros(period,24);
for T=1:period
    for t=1:24
        A(T,t)=vec(24*(T-1)+t);
    end
end

A_bar=mean(A,1);
A_bar_add=A_bar;

if period>1
    for T=2:period
        A_bar=[A_bar A_bar_add];
    end
end
figure
plot(A_bar)
hold on
plot(vec,'r')
title('total data: average tweets over time')
xlabel('time (hours)')
ylabel('# of tweets')
legend('average data','total data','location','NorthEastOutside')

% plot the difference of the total and the average data over time
Y_c=vec-A_bar;
figure
plot(Y_c)
title('total data: total minus average data over time')
xlabel('time (hours)')
ylabel('total - average')

% plot the difference at time t vs the difference at time t+1
figure
Y_c_next=Y_c;
Y=Y_c;
Y_c_next(1)=[];
Y(time)=[];
plot(Y,Y_c_next,'o')
title('total data: difference of total and average data at time t vs differnce at time t+1')
xlabel('difference at time t')
ylabel('difference at time t+1')

% make the same plots for each group
for g=1:G
   
   % split total data into one group
   initial=((g-1)*group_size)+1;
   final=g*group_size;
   group_mat=tweets_mat(initial:final,:);
   group_vec=sum(group_mat,1);
   
   % plot the total group data over time
   figure
   plot(group_vec)
   formatSpec= 'group %d: # of tweets over time';
   str=sprintf(formatSpec,g);
   title(str)
   xlabel('time (hours)')
   ylabel('# of tweets')
   
   % plot the total group data at time t vs the group data at time t+1
   figure
   A_c=group_vec;
   A_c_next=group_vec;
   A_c(time)=[];
   A_c_next(1)=[];
   plot(A_c,A_c_next,'o')
   formatSpec= 'group %d: # of tweets at time t vs # of tweets at time t+1';
   str=sprintf(formatSpec,g);
   title(str)
   xlabel('# of tweets at time t')
   ylabel('# of tweets at time t+1')
   
   % plot the average data for each group
   period=time/24;
   A=zeros(period,24);
   for T=1:period
       for t=1:24
           A(T,t)=group_vec(24*(T-1)+t);
       end
   end

   A_bar=mean(A,1);
   A_bar_add=A_bar;

   if period>1
       for T=2:period
           A_bar=[A_bar A_bar_add];
       end
   end
   figure
   plot(A_bar)
   hold on
   plot(group_vec,'r')
   formatSpec= 'group %d: average # of tweets over time';
   str=sprintf(formatSpec,g);
   title(str)
   xlabel('time (hours)')
   ylabel('# of tweets')
   legend('average data','total data','location','NorthEastOutside')
   
   % plot the difference of the groups total data and the groups average
   % data over time
   figure
   Y_c=group_vec-A_bar;
   plot(Y_c)
   formatSpec= 'group %d: total group data minus the average over time';
   str=sprintf(formatSpec,g);
   title(str)
   xlabel('time (hours)')
   ylabel('total - average')
   
   % plot the difference at time t vs the difference at time t+1
   figure
   Y_c_next=Y_c;
   Y=Y_c;
   Y_c_next(1)=[];
   Y(time)=[];
   plot(Y,Y_c_next,'o')
   formatSpec= 'group %d: difference of total and average data at time t vs difference at time t+1';
   str=sprintf(formatSpec,g);
   title(str)
   xlabel('difference at time t')
   ylabel('difference at time t+1')
end
end

function make_plots(tweets_mat,indep_vec,total_time,G,group_size,total_kicks_mat,indep_tweets_mat)
% plot of total data over time
%figure
% plot the tweets over time in black
vec=sum(tweets_mat,1);
%plot(vec,'k')
title('total data: # of tweets over time')
xlabel('time (hours)')
ylabel('number of tweets')
hold on
% plot the independent tweets over time for comparison in dashed line
%plot(indep_vec,'k--')

% make the same plot for individual groups
for g=1:2
   
   % split total data into one group
   initial=((g-1)*group_size)+1;
   final=g*group_size;
   group_mat=tweets_mat(initial:final,:);
   group_vec=sum(group_mat,1);
   indep_group_mat=indep_tweets_mat(initial:final,:);
   indep_group_vec=sum(indep_group_mat,1);
   
   % plot the total group data over time
   if g==1
       color='b';
   elseif g==2
       color='g';
   elseif g==3
       color='r';
   elseif g==4
       color='c';
   elseif g==5
       color='m';
   elseif g==6
       color='y';
   end
   
   str=sprintf('%s',color);
   str2=sprintf('%s--',color);
   plot(group_vec,str)
   plot(indep_group_vec,str2)
   legend('group 1','group 1 indep','kicks','group 2','group 2 indep','location','NorthEastOutside')
   
   y=[ylim];
   yl=y(2)
   line=total_kicks_mat(g,:)*yl;
   bar(line,str,'EdgeColor',str,'lineWidth',.00001,'barWidth',.005)
end
end

function make_plots_minus_seasonality(tweets_mat,indep_vec,total_time,G,group_size,total_kicks_mat)
% plot the tweets over time in black
vec=sum(tweets_mat,1);
hold on

% get the average data
period=total_time/24;
A=zeros(period,24);
for T=1:period
    for t=1:24
        A(T,t)=vec(24*(T-1)+t);
    end
end

A_bar=mean(A,1);
A_bar_add=A_bar;

if period>1
    for T=2:period
        A_bar=[A_bar A_bar_add];
    end
end

% plot the difference of the total and the average data over time in black
Y_c=vec-A_bar;
%figure
hold on
%plot(Y_c,'k')
title('minus seasonality')
xlabel('time (hours)')
ylabel('total - average')

% make the same plot for individual groups
for g=1:2
   
   % split total data into one group
   initial=((g-1)*group_size)+1;
   final=g*group_size;
   group_mat=tweets_mat(initial:final,:);
   group_vec=sum(group_mat,1);
   
   % plot the total group data over time
   if g==1
       color='b';
   elseif g==2
       color='g';
   elseif g==3
       color='r';
   elseif g==4
       color='c';
   elseif g==5
       color='m';
   elseif g==6
       color='y';
   end
   
   % get the average data for each group
   period=total_time/24;
   A=zeros(period,24);
   for T=1:period
       for t=1:24
           A(T,t)=group_vec(24*(T-1)+t);
       end
   end

   A_bar=mean(A,1);
   A_bar_add=A_bar;

   if period>1
       for T=2:period
           A_bar=[A_bar A_bar_add];
       end
   end
   
   % plot the difference of the groups total data and the groups average
   % data over time
   Y_c_group=group_vec-A_bar;
   
   str=sprintf('%s',color);
   plot(Y_c_group,str)
   legend('group 1','kicks','group 2','group 3','group 4','group 5','location','NorthEastOutside')
   
end 
   y=[ylim];
   yl=y(2)
   line=sum(total_kicks_mat)*yl;
   bar(line,.005)

end

function make_plots_random_groups(tweets_mat,indep_vec,total_time,G,group_size,N,total_kicks_mat)
% plot of total data over time
%figure
% plot the total tweets over time in black
vec=sum(tweets_mat,1);
%plot(vec,'k')
title('random groups: # of tweets over time')
xlabel('time (hours)')
ylabel('number of tweets')
hold on

for g=1:2
   % split total data into one group
   initial=((g-1)*group_size)+1;
   final=g*group_size;
   group_mat=tweets_mat(initial:final,:);
   group_vec=sum(group_mat,1);
   
   % plot the total group data over time
   if g==1
       plot(group_vec,'b')
   elseif g==2
       plot(group_vec,'r')
   end

end

% make random groups
rand_tweets_mat=tweets_mat(randsample(1:length(tweets_mat),length(tweets_mat)),:);
for w=1:2
   % split total data into one group
   initial=((w-1)*group_size)+1;
   final=w*group_size;
   rand_group_mat=rand_tweets_mat(initial:final,:);
   rand_group_vec=sum(rand_group_mat,1);
   
   % plot the random group data over time for 2 groups
   if w==1
       plot(rand_group_vec,'g')
   elseif w==2
       plot(rand_group_vec,'m')
   end
   legend('group 1','group 2','random group 1','random group 2','location','NorthEastOutside')
end
for t=1:total_time
    total_kicks_vec=sum(total_kicks_mat);
    if total_kicks_vec(t)>1
        total_kicks_vec(t)=1;
    end
end

   y=[ylim];
   yl=y(2)
   line=(total_kicks_vec)*yl;
   bar(line,.005)

end