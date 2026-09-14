include("input.jl")
using Base.Filesystem
using ITensors, ITensorMPS
using Base.Threads

let

  #===========================================================================#
             #================= Parameters ===================#
  #===========================================================================#
  Kitayev.N = parse(Int64, ARGS[1])
  Kitayev.theta = parse(Float64, ARGS[2])
  Kitayev.hz = parse(Float64, ARGS[3])
  # order = parse(Int64, ARGS[4])
  #Kitayev.hz = 1
  Kitayev.longrange = 1
  Kitayev.sites = siteinds("tJ",N;conserve_nf=true)
  #===========================================================================#
  #===========================================================================#

  Kitayev.J1=cos(pi*Kitayev.theta/2000)
  Kitayev.J2=sin(pi*Kitayev.theta/2000)

  Kitayev.t1=cos(pi*Kitayev.theta/2000)
  Kitayev.t2=sin(pi*Kitayev.theta/2000)
  #===========================================================================#
        #================= Printing of Parameters ===================#
  #===========================================================================#
  
  Kitayev.Values()
  Kitayev.printMatrix(Kitayev.JJ1)
  Kitayev.printMatrix(Kitayev.tt)

  #===========================================================================#
  H = Kitayev.Hamiltonian()  # Hamiltonian

  #===========================================================================#
     #================= DMRG to claculate Groundstate ===================#
  #===========================================================================#
  state = [isodd(n) ? "Up" : "Dn" for n=1:N]
  # state[4]="Emp"
  # state[15]="Emp"
  psi0 = MPS(Kitayev.sites,state)#random_mps(Kitayev.sites; linkdims=100)#
  @show flux(psi0)
  nsweeps = 20
  noise = [1e-3, 1e-6, 1e-8, 1e-10, 1e-12]
  maxdim = [10, 20, 100, 500, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]
  cutoff = [1E-10]
  energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
  println("\n Groundstate Energy: $energy \n")
  #===========================================================================#
  #===========================================================================#
  @show flux(psi)
  m000=0
  magz = expect(psi,"Sz")
  for (j,mz) in enumerate(magz)
    m000+=mz/N
      println("$j $mz")
  end
  println("Total Sz= $m000")
  #===========================================================================#
  #===========================================================================#
  file_name = "New"
  Kitayev.Make_Directory("../Data")
  Kitayev.Make_Directory("../Data/NSC")
  Kitayev.Make_Directory("../Data/NSC/MPS")
  # if order == 0
  #   directory_path = "../Data/NSC/SZ"
  #   Kitayev.Make_Directory(directory_path)
  # elseif order == 1
  #   directory_path = "../Data/NSC/SZiSZj"
  #   Kitayev.Make_Directory(directory_path)
  # end
  #===========================================================================#
  #===========================================================================#

  #===========================================================================#
                        #============Save MPS============#
  #===========================================================================#
  theta = Int(Kitayev.theta)
  hz = Int(Kitayev.hz * 100)
  using HDF5
  file_nameMPS = "../Data/NSC/MPS/SZ_$(Kitayev.N)_$(theta)_$(hz).h5"
  f = h5open(file_nameMPS,"w")
  write(f,"psi",psi)
  close(f)
  #===========================================================================#
  #===========================================================================#

  return
end
