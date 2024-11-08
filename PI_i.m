function f = PI_i(i,nx,N)

f = zeros(nx,N*nx);
f(:,(i-1)*nx+1:i*nx) = eye(nx);

return