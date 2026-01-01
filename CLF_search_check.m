size(Emax)
prod(size(Emax))
[m,I] = min(Emax,[],"all")

indx1 = floor(I/(length(p0)*length(p1)*length(thp)*length(z0)*length(z1)))+1
I = I-(indx1-1)*(length(p0)*length(p1)*length(thp)*length(z0)*length(z1))

indx2 = floor(I/(length(p0)*length(p1)*length(thp)*length(z0)))+1
I = I-(indx2-1)*(length(p0)*length(p1)*length(thp)*length(z0))

indx3 = floor(I/(length(p0)*length(p1)*length(thp)))+1
I = I-(indx3-1)*(length(p0)*length(p1)*length(thp))

indx4 = floor(I/(length(p0)*length(p1)))+1
I = I-(indx4-1)*(length(p0)*length(p1))

indx5 = floor(I/(length(p0)))+1
I = I-(indx5-1)*(length(p0))

indx6 = I


Emax(indx6,indx5,indx4,indx3,indx2,indx1)

