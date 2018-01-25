s = daq.createSession('ni');
s.IsContinuous = 1;
s.DurationInSeconds = 10;

s.Rate = 20000;
addAnalogInputChannel(s,'dev1', 'ai0', 'Voltage');
addAnalogInputChannel(s,'dev1', 'ai1', 'Voltage');

lh = addlistener(s,'DataAvailable', @plotMotorSpeed);


startForeground(s)
startBackground(s)


stop(s)

[data,timeStamps,triggerTime] = startForeground(s1);
a=[1,2;3,4]
msg = [a,b;['yes'] ['no']]
a=123;b=456;
disp([num2str(a), num2str(b) ,'yes';1,2,3])
	
X = rand(1,3);
disp('Corn	Oats	Hay'\n)
disp(X)