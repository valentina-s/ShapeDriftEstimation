function nothing = test2()
temp = zeros(100,1000,10);
result = zeros(100,1);
parfor k=1:100
    result(k,1) = fun(temp(k,:,:));
end
nothing = 1
end

function result = fun(temp)
    result = sum(sum(temp))
end
    



