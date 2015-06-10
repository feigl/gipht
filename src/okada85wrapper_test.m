load fitfun28
kp = 9;
uokada1 = okada85_wrapper(pg(kp:kp+9),[DST.x';DST.y'],nu);
quiver(DST.x',DST.y',uokada1(1,:),uokada1(2,:));