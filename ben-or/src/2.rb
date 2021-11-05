def seki(a , b)

    c=0;
    while(a!=0)
    if ((a & 1)==1)
    c^=b;
    end
    b<<=1; a>>=1;
    end
    
return c;

end

a=111111111111111111111111111111111
b=2777777777777777777777777777777779
print(a.to_s(2))
print("\n")
print(b.to_s(2))
print("\n")
print(seki(a,b).to_s(2))
print("\n")
