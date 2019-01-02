function df_all = df_all(df,x,n,N)

df_all = zeros((N+1)*n);
for j = 1:N+1
    df_all(1+(j-1)*n:j*n,1+(j-1)*n:j*n) = df(x(1+(j-1)*n:j*n));
end
    
    

end