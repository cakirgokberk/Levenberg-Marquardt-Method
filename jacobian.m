function J = jacobian(a,b,x)

syms xi ai bi

J_var = jacobian(log10(ai*xi)+3*cos(bi*xi),[ai, bi]);
xi = transpose(x);
ai = a;
bi = b;
J = double(subs(J_var));

end

