%this script runs 1000 simulations for planform meander river development
%with a range of channel bed slope
R = 1.65;%submerged specific gravity of sand
v = 1e-6;%viscosity
g = 9.81;%gravitational acceleration

noSimulation  = 1000;
channelSlopeL = log10(6e-5);%upper limit for slope input
channelSlopeU = log10(1e-3);%lower limit for slope input
sInput = 10.^(channelSlopeL + (channelSlopeU-channelSlopeL)...
    * rand(1,noSimulation));
sInput = sort(sInput);

CfInput = 0.1975 * sInput.^0.4068;%friction coefficient
hInput = 1220 .* sInput.^(-0.47) * (R*v)^(2/3)/g^(1/3);%flow depth
bInput = 20.23 .* hInput.^1.266/2;%half channel width
uInput = 35 .* sInput.^(0.26) * (R*g*v)^(1/3) ./ sqrt(CfInput);%reach-averaged flow velocity
EInput = fliplr((1:9/noSimulation:10)*2e-7);%erosion coefficient

% parfor i =1:noSimulation%paralell computing
for i =1:noSimulation %use for loop if parfor is not available
    rTmp = simuFunction(hInput(i),bInput(i),uInput(i),CfInput(i),EInput(i),sInput(i)); 
    Result(i).r = rTmp;
    fprintf('No. of simulation is %d \n' ,i);
end





