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
  #order = parse(Int64, ARGS[7])
  Kitayev.hz = parse(Float64, ARGS[6])
  Kitayev.longrange = 1
  Kitayev.sites = siteinds("tJ",N;conserve_nf=true)
  #===========================================================================#
  #===========================================================================#

  #===========================================================================#
        #================= Printing of Parameters ===================#
  #===========================================================================#
  
  Kitayev.Values()
  # Kitayev.printMatrix(Kitayev.JJ1)
  # Kitayev.printMatrix(Kitayev.tt)

  #===========================================================================#
  H = Kitayev.Hamiltonian()  # Hamiltonian

  #===========================================================================#
     #================= DMRG to claculate Groundstate ===================#
  #===========================================================================#
  state = [isodd(n) ? "Up" : "Dn" for n=1:N]
  # state[2]="Emp"
  # state[9]="Emp"
  psi0 = MPS(Kitayev.sites,state)#random_mps(Kitayev.sites; linkdims=100)#
  @show flux(psi0)
  nsweeps = 20
  noise = [1e-3, 1e-6, 1e-8, 1e-10, 1e-12]
  maxdim = [10, 20, 100, 500, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]
  cutoff = [1E-10]
  energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)#, noise
  println("\n Groundstate Energy: $energy \n")
  #===========================================================================#
  #===========================================================================#

   #===========================================================================#
  #===========================================================================#
  file_name = "New"
  Kitayev.Make_Directory("../Data")
  Kitayev.Make_Directory("../Data/NSC")
  Kitayev.Make_Directory("../Data/NSC/Test/")
  #===========================================================================#
  #===========================================================================#

  #println("$ii \t $jj")
  psi00 = psi
  psi01 = psi
  psi1 =psi#apply(O1,psi;cutoff=1e-12,maxdim=1000) # O_ii|g⟩= |ii⟩
  psi2 =psi#apply(O2,psi00;cutoff=1e-12,maxdim=1000) # O_jj|g⟩= |jj⟩
  print(inner(psi2,psi1), "\n")

  #===========================================================================#
  #===========================================================================#
  J1 = Int(Kitayev.J1 * 100)
  J2 = Int(Kitayev.J2 * 100)
  t1 = Int(Kitayev.t1 * 100)
  t2 = Int(Kitayev.t2 * 100)
  hz = Int(Kitayev.hz * 100)
  
  file_name = "../Data/NSC/Test/SZ_$(Kitayev.N)_$(J1)_$(J2)_$(t1)_$(t2)_$(hz).dat"
  #===========================================================================#
  #===========================================================================#
  file = open(file_name, "w")
  println(file, "#ttotal <Sz(SL,ttotal).Sz(SR,0)> t(ms) SR   SL")
  #===========================================================================#
           #================= Time Evolution ===================#
  #===========================================================================#

  #===========================================================================#
  #===========================================================================#
    cutofft = 1E-8
    tau = 0.01
    ttotal = 30.0
    gates=Kitayev.TEBD_gates(tau)

    for t in 0.0:tau:ttotal
      x = psi1 # O_ii|g⟩= |ii⟩
      y = psi2 # O_jj|g⟩= |jj⟩
      prefactor = exp(im * 0 * energy) # exp(+iE_gt)
      Sz = inner(x',H,x)
      Sz = Sz * prefactor # exp(+iE_gt)x⟨g| O_ii exp(-iHt) O_jj |g⟩
      SZreal = real(Sz-energy)
      SZimag = imag(Sz)
      println(file, "$t \t $SZreal \t $SZimag")
      println("$t \t $SZreal \t $SZimag")
      t ≈ ttotal && break
      psi1 = apply(gates, psi1; cutoff=1e-12, maxdim=1000) # Gate Evolution exp(-iHt)|jj⟩
    end
#===========================================================================#
#===========================================================================#

 

  return
end