function testMexTriDTriD

[Qt,c,sDim,Js] = mexDiagTriDQOP(3,3,2000);

full(Qt)

full(c)

sDim

Js

return