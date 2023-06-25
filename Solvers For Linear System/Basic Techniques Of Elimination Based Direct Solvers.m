classdef functionscontainer
    methods

        %Function to solve using forward substitution
        function [u] = FdSubs(obj,A,a)
        [n,n] = size(A);
        u = zeros(n,1);
        for t=1:n
            i=n+1-t;
            u(i)=(a(i)-(A(i,:)*u))/(A(i,i));
        end
        end
        
        %Function to solve using backward substitution
        function [v] = myans(obj,A,a)
        [n,n] = size(A);
        v = zeros(n,1);
        for t=1:n
            j=n+1-t;
               v(j)=(a(j)-(A(j,:)*v))/(A(j,j));
        end
        end

        %Function for LU decomposition
        function [L,U] = mylu(obj,A)
        n = size(A,1);
        for k = 1:n
            if A(k,k)==0
                warning('LU factorization cannot be done');
                return; 
            end
            i = k+1:n;
            A(i,k) = A(i,k)/A(k,k);
            A(i,i) = A(i,i)-A(i,k)*A(k,i);
        end
        L = tril(A, -1) + eye(n);
        U = triu(A);
        end
        
        % LU factorization with partial (row) pivoting
        function [L,U,P] = PLU_pivot(obj,A)
        [n,n]=size(A);
        L=eye(n); P=L; U=A;
        for k=1:n
            [pivot m]=max(abs(U(k:n,k)));
            m=m+k-1;
            if m~=k
                % interchange rows m and k in U
                temp=U(k,:);
                U(k,:)=U(m,:);
                   U(m,:)=temp;
                % interchange rows m and k in P
                temp=P(k,:);
                P(k,:)=P(m,:);
                P(m,:)=temp;
                if k >= 2
                    temp=L(k,1:k-1);
                    L(k,1:k-1)=L(m,1:k-1);
                    L(m,1:k-1)=temp;
                end
            end
            for j=k+1:n
                L(j,k)=U(j,k)/U(k,k);
                U(j,:)=U(j,:)-L(j,k)*U(k,:);
            end
        end
        end
    end
end