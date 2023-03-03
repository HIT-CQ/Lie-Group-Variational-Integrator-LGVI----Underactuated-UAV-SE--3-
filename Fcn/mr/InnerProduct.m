function IP = InnerProduct(A,B) 
% Input:  A (n×m向量)        B (n×m向量)       
 % 
% Output: <A,B> 


IP = trace(A.'*B);
    
end