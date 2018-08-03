function vel=circles(x,y),

u=4*x-2; v=4*y-2;

vel=tanh(30*(u.^2+v.^2-0.25)) ...
	+tanh(30*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30*((u+0.75).^2+(v+0.75).^2-0.25)) ;
