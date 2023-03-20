function testArrowTriDQOP
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sDim = 3; kBlock = 3; randSeed = 2009;

[Qt,c, sDim, Js] = mexArrowTriDQOP(sDim,kBlock,randSeed);

full(Qt)

% full(c)

full(sDim)

full(Js)

end

