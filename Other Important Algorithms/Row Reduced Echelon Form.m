n = 3;
A = rand(n)
%A(1,1) = 100;
X = row_echelon(A)
Y = rref(A)

function [row_echelon_matrix] = row_echelon(A) 
    [m, n] = size(A);  % get the size of the matrix 
    row_echelon_matrix = A;  % initialize the row echelon matrix as the original matrix 
     
    for i = 1:m  % loop through the rows 
        for j = 1:n  % loop through the columns 
            if row_echelon_matrix(i,j) ~= 0  % if the element is non-zero 
                % divide the row by the element to make the leading coefficient 1 
                row_echelon_matrix(i,:) = row_echelon_matrix(i,:) / row_echelon_matrix(i,j); 
                break;  % move on to the next row 
            end 
        end 
         
        for k = 1:m  % loop through the rows again 
            if k ~= i  % skip the current row 
                % subtract a multiple of the current row from the other row to make the leading coefficient 0 
                row_echelon_matrix(k,:) = row_echelon_matrix(k,:) - row_echelon_matrix(k,j) * row_echelon_matrix(i,:); 
            end 
        end 
    end 
end 