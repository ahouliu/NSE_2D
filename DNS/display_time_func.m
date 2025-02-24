function [time] = display_time_func()

time = strsplit(string(datetime));
time = time(2);
fprintf('Now is %s.\n',time);

end