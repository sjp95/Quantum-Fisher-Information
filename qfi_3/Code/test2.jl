include("input.jl")
using Base.Filesystem
using ITensors, ITensorMPS
using Base.Threads

let

  #===========================================================================#
             #================= Parameters ===================#
  #===========================================================================#
  Kitayev.N = parse(Int64, ARGS[1])
  Kitayev.J1 = parse(Float64, ARGS[2])
  Kitayev.J2 = parse(Float64, ARGS[3])
  Kitayev.t1 = parse(Float64, ARGS[4])
  Kitayev.t2 = parse(Float64, ARGS[5])
  order = parse(Int64, ARGS[6])
  ii=parse(Int64, ARGS[7])
  jj=parse(Int64, ARGS[8])
  Kitayev.hz = parse(Float64, ARGS[9])
  Kitayev.longrange = 1
  #===========================================================================#
  #===========================================================================#

  #===========================================================================#
        #================= Printing of Parameters ===================#
  #===========================================================================#
  
  Kitayev.Values()
  # Kitayev.printMatrix(Kitayev.JJ1)
  # Kitayev.printMatrix(Kitayev.tt)

  #===========================================================================#
  #===========================================================================#
                    #============Reading MPS============#
  #===========================================================================#
  J1 = Int(Kitayev.J1 * 100)
  J2 = Int(Kitayev.J2 * 100)
  t1 = Int(Kitayev.t1 * 100)
  t2 = Int(Kitayev.t2 * 100)
  hz = Int(Kitayev.hz * 100)

  using HDF5
  file_nameMPS = "../Data/NSC/MPS/SZ_$(Kitayev.N)_$(J1)_$(J2)_$(t1)_$(t2)_$(hz).h5"
  f = h5open(file_nameMPS,"r")
  psi = read(f,"psi",MPS)
  close(f)
  #===========================================================================#
  #===========================================================================#
  Kitayev.sites = siteinds(psi)
  H = Kitayev.Hamiltonian()  # Hamiltonian
  #===========================================================================#
     #================= DMRG to claculate Groundstate ===================#
  #===========================================================================#
  @show flux(psi)
  
  energy=inner(psi',H,psi)
  println("\n Groundstate Energy: $energy \n")
  #===========================================================================#
  #===========================================================================#

  #===========================================================================#
  #===========================================================================#
  file_name = "New"
  Kitayev.Make_Directory("../Data")
  Kitayev.Make_Directory("../Data/NSC")
  Kitayev.Make_Directory("../Data/NSC/MPS")
  if order == 0
    directory_path = "../Data/NSC/SZ"
    Kitayev.Make_Directory(directory_path)
  elseif order == 1
    directory_path = "../Data/NSC/SZiSZj"
    Kitayev.Make_Directory(directory_path)
  elseif order == 2
    directory_path = "../Data/NSC/SZiSZjSZk"
    Kitayev.Make_Directory(directory_path)
  end
  #===========================================================================#
  #===========================================================================#
    O1 = Kitayev.Operators(order, ii)
    O2 = Kitayev.Operators(order, jj)

    println("$ii \t $jj ")
    psi00 = psi
    psi01 = psi
    psi1 = apply(O1, psi01; cutoff=1e-12, maxdim=1000) # O_ii|g⟩= |ii⟩
    psi2 = apply(O2, psi00; cutoff=1e-12, maxdim=1000) # O_jj|g⟩= |jj⟩
    #println(inner(psi2, psi1))

    

    #===========================================================================#
    #===========================================================================#
    if order == 0
      file_name = "../Data/NSC/SZ/SZ_$(Kitayev.N)_$(J1)_$(J2)_$(t1)_$(t2)_$(ii)_$(jj)_$(hz).dat"
    elseif order == 1
      file_name = "../Data/NSC/SZiSZj/SZiSZj_$(Kitayev.N)_$(J1)_$(J2)_$(t1)_$(t2)_$(ii)_$(jj)_$(hz).dat"
    elseif order == 2
      file_name = "../Data/NSC/SZiSZjSZk/SZiSZjSZk_$(Kitayev.N)_$(J1)_$(J2)_$(t1)_$(t2)_$(ii)_$(jj)_$(hz).dat"
    end
    #===========================================================================#
    #===========================================================================#
    file = open(file_name, "w")
    println(file, "#ttotal <Sz(SL,ttotal).Sz(SR,0)> t(ms) SR   SL")
    #===========================================================================#
    #===========================================================================#

    cutofft = 1E-8
    tau = 0.01
    ttotal = 50.0
    gates=Kitayev.TEBD_gates(tau)

    for t in 0.0:tau*10:ttotal
      x = psi1 # O_ii|g⟩= |ii⟩
      y = psi2 # O_jj|g⟩= |jj⟩
      prefactor = exp(im * t * energy) # exp(+iE_gt)
      Sz = inner(y, x)
      Sz = Sz * prefactor # exp(+iE_gt)x⟨g| O_ii exp(-iHt) O_jj |g⟩
      SZreal = real(Sz)
      SZimag = imag(Sz)
      println(file, "$t \t $SZreal \t $SZimag  \t $ii \t $jj ")
      println("$t \t $SZreal \t $SZimag")
      t ≈ ttotal && break
      for i in 1:10
        psi1 = apply(gates, psi1; cutoff=1e-12, maxdim=1000) # Gate Evolution exp(-iHt)|jj⟩
      end
    end

  return
end
