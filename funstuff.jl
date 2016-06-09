N=10^6
v=linspace(0.1,100/3,N)

function forloop(N::Int,v::Vector)
  for i=1:N
    a=v[i]
    b=a<=100/3
  end
  end

function forloop1(N::Int,v::Vector)
  for i=1:N
    a=v[i]
  end
  end


function whileloop(N)
    v=1.1
  while  v<= 100/3
    a=v
    v+=(100/3-1.1)/N
    end
  end


function whileloop1(N::Int,v::Vector)
    i=1;
  while  v[i]< 100/3
     a=v[i]
    i+=1;
    end
  end

@time forloop(N,v)
@time forloop1(N,v)
@time whileloop(N)
@time whileloop1(N,v)


#global vars

a=2.5;
b=9;

function globvars(s)
  y=s*a+b;
end


function locvars(s,a,b)
  y=s*a+b;
end

function taketimeg(N)
  for i=1:N
    b=globvars(i)
  end
end

function taketimel(N)
   for i=1:N
    c=locvars(i,a,b)
  end
end

function taketimef(N)
   for i=1:N
    function func(s)
      y=s*a+b;
    end
    c=locvars(i,a,b)
  end
end

@time taketimeg(10^7)
#1.4 seconds
@time taketimel(10^7)
#1.4 seconds
@time taketimef(10^7)
#3.1 seconds
