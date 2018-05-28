function [G_red,T]=CLLL(H)
% complex LLL algorithm
[Q,R] = qr(H);
K = size(H,2);
T = eye(K);
delta = 0.75;  % (1/2,1)
k = 2;
while k <= K
    for n = k-1:-1:1
        mu = round(R(n,k)/R(n,n));
        if mu~=0
            R(1:n,k) = R(1:n,k)-mu*R(1:n,n);
            T(:,k) = T(:,k)-mu*T(:,n);
        end
    end
    if (delta*norm(R(k-1,k-1))^2 > norm(R(k,k))^2+norm(R(k-1,k))^2)  %列变换
        R(:,[k-1 k]) = R(:,[k k-1]);
        T(:,[k-1 k]) = T(:,[k k-1]);
        %计算Givens旋转矩阵，以使R(k,k-1)为0
%         alpha = R(k-1,k-1)/sqrt(R(k-1:k,k-1).'*R(k-1:k,k-1));
%         beta = R(k,k-1)/sqrt(R(k-1:k,k-1).'*R(k-1:k,k-1));
        alpha = R(k-1,k-1)/norm(R(k-1:k,k-1));
        beta = R(k,k-1)/norm(R(k-1:k,k-1));
        Theta = [alpha' beta';-beta alpha];
        R(k-1:k,k-1:K) = Theta*R(k-1:k,k-1:K);
        Q(:,k-1:k) = Q(:,k-1:k)*Theta';
        k = max([k-1 2]);
    else
        k = k+1;
    end
end
G_red = Q*R;



