function [x]=log2matSHORT(file)

[trial, event_type, code, time] = textread(file, '%*s %d %s %s %f %*[^\n]', 'headerlines', 5);%---Liest das Logfile ein;
x{:,:}=trial;
x{:,2}=event_type;
x{:,3}=code;
x{:,4}=time;
