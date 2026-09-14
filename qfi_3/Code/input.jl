#using ITensors, ITensorMPS
module Kitayev
using Base.Filesystem
using ITensors, ITensorMPS
#===========================================================================#
             #================== Parameters ==================#
#===========================================================================#
export N, J1, J2, t1, t2, JJ1, tt,hz, longrange, sites

N = 4  # Change this to the desired size
J1=1
J2=0.5
t1=1
t2=0.5
hz=500
theta=0
longrange=1
#sites = siteinds("tJ",N;conserve_qns=false)

JJ1 = zeros(Float64, N, N)
tt = zeros(Float64, N, N)

sites = siteinds("tJ",N;conserve_qns=false)
gates = ITensor[]

#===========================================================================#
    #================== Function to take the values ==================#
#===========================================================================#
function Resize()
        global JJ1 = zeros(Float64, N, N)
        global tt = zeros(Float64, N, N)
    end
#===========================================================================#
#===========================================================================#
function Values()
    Resize()
    for i in 1:3:N
        #===============#
        #===============#
        tt[i,i+1]=t1
        tt[i+1,i]=t1

        JJ1[i,i+1]=J1
        JJ1[i+1,i]=J1
        #===============#
        #===============#
        tt[i+1,i+2]=t1
        tt[i+2,i+1]=t1

        JJ1[i+1,i+2]=J1
        JJ1[i+2,i+1]=J1
        #===============#
        #===============#
        if i+3<=N
            tt[i+2,(i+2)%N+1]=t2
            tt[(i+2)%N+1,i+2]=t2

            JJ1[i+2,(i+2)%N+1]=J2
            JJ1[(i+2)%N+1,i+2]=J2
        end
           
    end

end
#===========================================================================#
#===========================================================================#
function printMatrix(H0)
    rows, cols = size(H0)
    for i in 1:rows
        for j in 1:cols
          print(H0[i,j], "\t")
        end
        print("\n")
      end
print("\n")
end
#===========================================================================#
#===========================================================================#

function Make_Directory(directory_path)
    if !isdir(directory_path)
        # Create the directory
        mkdir(directory_path)
    else
        println("Directory already exists: $directory_path")
    end
end
#===========================================================================#
#===========================================================================#

# #===========================================================================#
#       #================== Hamiltonian as an MPO ==================#
# #===========================================================================#
function Hamiltonian()
    #global sites = siteinds("tJ",N;conserve_qns=false)
    os = OpSum()
    for j=1:N-1
        os += JJ1[j,j+1],"Sz",j,"Sz",j+1
        os += JJ1[j,j+1]/2.0,"S+",j,"S-",j+1
        os += JJ1[j,j+1]/2.0,"S-",j,"S+",j+1

        os += -JJ1[j,j+1]/4.0,"Ntot",j,"Ntot",j+1 #// Hopping term

        os += -tt[j,j+1],"Cdagup",j,"Cup",j+1 #// Hopping term
        os += -tt[j,j+1],"Cdagup",j+1,"Cup",j
        os += -tt[j,j+1],"Cdagdn",j,"Cdn",j+1
        os += -tt[j,j+1],"Cdagdn",j+1,"Cdn",j
        os += -hz*1.0,"Sz",j,"I",j+1
    end
    println("hz=$hz")
    os += -hz*1.0,"Sz",N,"I",N-1
    H = MPO(os,sites)
    return H
end
#===========================================================================#
#===========================================================================#


#====================================================================================#
          #================== Operators as an MPO ==================#
#====================================================================================#
function Operators(order,ii)
    os1 = OpSum()

    for i in ii-longrange:ii+longrange



        if order==0


            #---------------------------#
            if i==ii-longrange
                os1 += 1.0,"Sz", ii
            end
            #---------------------------#


        elseif order ==1

            #-------------------------------------------------------------------
            if i >= 1 && i <= N
                if abs(JJ1[ii, i]) > 0.0
                #---------------------------------------------------------------
                os1 += JJ1[ii, i]*0.5,"Sx",ii,"S+",ii,"S+",i; 
                os1 += JJ1[ii, i]*0.5,"Sx",ii,"S-",ii,"S+",i;

                os1 += JJ1[ii, i],"Sx",ii,"Sz",ii,"Sz",i;
                #---------------------------------------------------------------
                end
            end
            #------------------------------------------------------------------

        elseif order ==2


            for kpp in ii-longrange:ii+longrange

                if kpp >= 1 && kpp <= N
                    #-------------------------------------------------------------------
                    if i >= 1 && i <= N
                    
                        if abs(JJ1[ii, i]*JJ1[ii,kpp]) > 0.0
                            #--------------------------------------------------------------------------------
                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.25,"Sx",ii,"S+",ii,"S-",i,"S+",ii,"S-",kpp; 
                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.25,"Sx",ii,"S+",ii,"S-",i,"S-",ii,"S+",kpp; 
                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.5,"Sx",ii,"S+",ii,"S-",i,"Sz",ii,"Sz",kpp; 

                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.25,"Sx",ii,"S-",ii,"S+",i,"S+",ii,"S-",kpp; 
                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.25,"Sx",ii,"S-",ii,"S+",i,"S-",ii,"S+",kpp; 
                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.5,"Sx",ii,"S-",ii,"S+",i,"Sz",ii,"Sz",kpp; 

                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.5,"Sx",ii,"Sz",ii,"Sz",i,"S+",ii,"S-",kpp;
                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*0.5,"Sx",ii,"Sz",ii,"Sz",i,"S-",ii,"S+",kpp;
                            os1 += JJ1[ii, i]*JJ1[ii,kpp]*1.0,"Sx",ii,"Sz",ii,"Sz",i,"Sz",ii,"Sz",kpp;
                            #--------------------------------------------------------------------------------
                        end

                    end
                    #-------------------------------------------------------------------
                end

            end



        end


    end

    O1=MPO(os1,sites)
    return O1
end
#====================================================================================#
#====================================================================================#




#===========================================================================#
            #================== TEBD GATES ==================#
#===========================================================================#
function TEBD_gates(tau)

    gates = ITensor[]

    for j in 1:(N - 1)
        s1 = sites[j]
        s2 = sites[j + 1]

        hj = JJ1[j,j+1] * op("Sz", s1) * op("Sz", s2) +
             JJ1[j,j+1]/2.0 * op("S+", s1) * op("S-", s2) +
             JJ1[j,j+1]/ 2.0 * op("S-", s1) * op("S+", s2) -
             JJ1[j,j+1]/4.0 * op("Ntot",s1) * op("Ntot",s2) -
             tt[j,j+1] * op("Cdagup",s1)* op("Cup",s2) - 
             tt[j,j+1] * op("Cdagup",s2)* op("Cup",s1) -
             tt[j,j+1] * op("Cdagdn",s1)* op("Cdn",s2) -
             tt[j,j+1] * op("Cdagdn",s2)* op("Cdn",s1)- hz * op("Sz", s1)*op("Id", s2)

             if j==N-1
                 hj+=- hz * op("Sz", s2)*op("Id", s1)
             end

        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
      end
      #   s1 = sites[N-1]
      #   s2 = sites[N]
      #   hj =- hz * op("Sz", s2)*op("Id", s1)
      #   Gj = exp(-im * tau / 2 * hj)
      #   push!(gates, Gj)
      # Include gates in reverse order too
      # (N,N-1),(N-1,N-2),...
      append!(gates, reverse(gates))
      return gates
end
#===========================================================================#
#===========================================================================#




end
using .Kitayev
#using .Kitayev.Operators
