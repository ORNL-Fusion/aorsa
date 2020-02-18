% Read Amat, Brhs, and solution vector from files for Azzam to check

T = readtable('p_brhs_for_azzam_ed.txt','ReadVariableNames',false,'HeaderLines',1);
b = complex(table2array(T(:,2)),table2array(T(:,3)));

T = readtable('p_solution_for_azzam_ed.txt','ReadVariableNames',false,'HeaderLines',1);
x = complex(table2array(T(:,2)),table2array(T(:,3)));

T = readtable('p_mat_for_azzam_ed.txt','ReadVariableNames',false,'HeaderLines',1);
ai = table2array(T(:,1));
aj = table2array(T(:,2));
Aij = complex(table2array(T(:,3)),table2array(T(:,4)));
N = numel(ai);
for i=1:numel(ai)
    A(ai(i),aj(i)) = Aij(i);
end

% check norm(A*x-b)
disp(['N : ', num2str(numel(b))]);
disp(['norm_1(A) : ', num2str(norm(A,1))]);
disp(['norm_inf(A) : ', num2str(norm(A,Inf))]);

l2 = norm(A*x-b);
disp(['norm_2(A*x-b) : ', num2str(l2)]);

cn = cond(A);
disp(['cond_2(A) : ', num2str(cn)]);
disp(['rcond_1(A) : ', num2str(rcond(A))]);

err1 = l2 / (norm(A)*norm(b));
disp(['norm_2(Ax-b)/(norm_2(A)*norm_2(b)) : ', num2str(err1)]);

figure(1)
plot(real(A*x),real(b))
figure(2)
plot(imag(A*x),imag(b))


disp('done')