using LinearAlgebra
using ColorSchemes, LaTeXStrings, Colors
using CairoMakie
using ProgressMeter
using Printf
using HDF5

const NUMBLASTHREADS = 10
const COLOUR = :imola
const CMAP = cgrad(COLOUR)
const LSTYLES = [:solid :dash :dot :dashdotdot]

function qfi_bounds(N::Int)
  bounds = zeros(Float64, N + 1)
  for k in 1:N
    s = fld(N, k)               # floor(N/k)
    bounds[k+1] = s * k^2 + (N - s * k)^2
  end
  return bounds
end

function hm_Tr_cpp(n=2, thn=3, h::Float64=0.0;
                        exe          = "../build/trun",
                        out::Symbol  = :f,
                        verbose::Bool = true,
                        clabels=true)

    N        = 3n
    filename = "results/trimerQFI_hm_Tr$(n)_h$(round(h, digits=2)).png"
    theta_values = range(0, pi/2, length=2^thn+1)
    T_values     = exp.(range(log(1e-4), log(3.0), length=2^(thn+2)+1))
#     println(T_values)

    Tmin = T_values[1]
    Tmax = T_values[end]
    nT   = length(T_values)

    #------------------------------------------------------------------
    # Helper: match C++ filename convention exactly
    # round(x * 100) as integer, same as static_cast<int>(std::round(x*100))
    #------------------------------------------------------------------
    fmt(x) = lpad(round(Int, x * 1e6), 8, '0')

    #------------------------------------------------------------------
    # Julia-side HDF5 cache — keyed by (N, h, ntheta, nT)
    # Stores the full qfi_matrix so re-runs skip all computation
    #------------------------------------------------------------------
    cachefile = "results/cache_hm_Tr_N$(N)_h$(fmt(h))_nth$(length(theta_values))_nT$(nT).h5"
    mkpath("results")

    recomp = 1
    # Try loading from cache first
    if isfile(cachefile)
        global qfi_matrix, T_vals_cached, theta_vals_cached  # Declare upfront
        try
            verbose && println("Cache found — loading qfi_matrix from $cachefile")
            h5open(cachefile, "r") do f
                qfi_matrix = read(f, "qfi_matrix")
                T_vals_cached     = read(f, "T_values")
                theta_vals_cached = read(f, "theta_values")

                # Sanity check dimensions match current call
                if size(qfi_matrix) == (length(theta_values), length(T_values))
                    verbose && println("Cache loaded successfully — skipping C++ calls.")
                    recomp=0
#                     return (qfi, T_vals, theta_vals)
                else
                    @warn "Cache dimensions mismatch — recomputing."
                    recomp=1
                end
            end
#             return (qfi, T_vals, theta_vals)
        catch e
            @warn "Cache read failed: $e — recomputing."
            recomp=1
        end
    end
    if recomp == 0
        hmap = _plot_heatmap(qfi_matrix, T_vals_cached, theta_vals_cached,
                                        N, h, filename, out, clabels)
        return hmap
    end

    #------------------------------------------------------------------
    # Cache miss — run C++ and collect qfi_matrix
    #------------------------------------------------------------------
    qfi_matrix = zeros(length(theta_values), length(T_values))

    p = verbose ? Progress(length(theta_values); dt=0.5, barlen=50, color=:yellow) : nothing

    for (ir, theta) in enumerate(theta_values)

        J1 = cos(theta)
        J2 = sin(theta)
        hz = h

        #--------------------------------------------------------------
        # Call C++ executable once per theta — sweeps all T internally
        # Arguments: N  J1*100  J2*100  hz*100  Tmin  Tmax  nT
        #--------------------------------------------------------------
#         cmd = `$exe $N $(J1) $(J2) $(hz) $Tmin $Tmax $nT`
#         fmt_arg(x) = @sprintf("%.6f", x)
        fmt_arg(x) = round(x, digits=6)
        cmd = `$exe $N $(fmt_arg(J1)) $(fmt_arg(J2)) $(fmt_arg(hz)) $(fmt_arg(Tmin)) $(fmt_arg(Tmax)) $nT 1`
#         println(cmd)
        try
#             run(pipeline(cmd))
            run(pipeline(cmd, stdout=devnull, stderr=devnull))
        catch e
            @warn "C++ call failed for theta=$(round(theta, digits=4)): $e"
            verbose && next!(p)
            continue
        end

        #--------------------------------------------------------------
        # One dat file per (N, J1, J2, hz) — all T rows inside
        # Format per line: T hz J1 J2 qfi_density
        #--------------------------------------------------------------
        datfile = "../Data/QFI/QFI_$(N)_$(fmt(J1))_$(fmt(J2))_$(fmt(hz)).dat"

        if !isfile(datfile)
            @warn "Missing output file: $datfile"
            verbose && next!(p)
            continue
        end

        lines = readlines(datfile)
        if isempty(lines)
            @warn "Empty file: $datfile"
            verbose && next!(p)
            continue
        end

        # Parse all lines into (T => qfi) dict then match to T_values
        # by nearest neighbour — avoids any float formatting mismatch
        parsed = Dict{Float64,Float64}()
        for l in lines
            isempty(strip(l)) && continue
            cols = parse.(Float64, split(strip(l)))
            parsed[cols[1]] = cols[5]   # T => qfi_density
        end

        cpp_T_vals = collect(keys(parsed))
        for (iT, Tv) in enumerate(T_values)
            nearest = cpp_T_vals[argmin(abs.(cpp_T_vals .- Tv))]
            qfi_matrix[ir, iT] = parsed[nearest]
        end

        verbose && next!(p)
    end

    #------------------------------------------------------------------
    # Write cache — stores everything needed to reproduce the plot
    #------------------------------------------------------------------
    try
        h5open(cachefile, "w") do f
            write(f, "qfi_matrix",   qfi_matrix)
            write(f, "T_values",     collect(T_values))
            write(f, "theta_values", collect(theta_values))
        end
        verbose && println("Cache written to $cachefile")
    catch e
        @warn "Cache write failed: $e"
    end

    return _plot_heatmap(qfi_matrix, collect(T_values), collect(theta_values),
                         N, h, filename, out, clabels)
end


#----------------------------------------------------------------------
# Plotting extracted into a helper so both the cache-hit and
# cache-miss paths call identical code
#----------------------------------------------------------------------
function _plot_heatmap(qfi_matrix, T_values, theta_values, N, h, filename, out, clabels=false)

    tick_ratios = [0, 0.5, 1, 2, Inf]
    tick_pos    = atan.(tick_ratios)
    tick_labels = ["0", "0.5", "1", "2", "∞"]

    log_T_values = log10.(T_values)
    log_T_ticks  = floor(log_T_values[1]):1:ceil(log_T_values[end])
    mat          = Float64.(qfi_matrix)

    fig = Figure(size=(400, 300))
    ax  = Axis(fig[1, 1],
               xlabel = L"J_2/J_1",
               ylabel = L"$T$ (in log scale)",
               xticks = (tick_pos, tick_labels),
               yticks = (log_T_ticks, [L"10^{%$(Int(e))}" for e in log_T_ticks]))
    xlims!(ax, (theta_values[1], theta_values[end]))

    hm = heatmap!(ax,
                  collect(theta_values),
                  log_T_values,
                  mat,
                  colormap   = COLOUR,
                  colorrange = (0.0, 4.0)
                  )

    text!(ax, 0.02, 0.98,
          text  = L"h=%$(round(h, digits=2))",
          space = :relative, color = :white,
          align = (:left, :top))

    if clabels==true
        contour!(ax,
                collect(theta_values),
                log_T_values,
                mat,
                levels    = Float64.(qfi_bounds(N) ./ N),
                labels    = true,
                labelsize = 16,
                color     = :black,
                linestyle = :dash,
                linewidth = 1.5)
    else
        contour!(ax,
                collect(theta_values),
                log_T_values,
                mat,
                levels    = Float64.(qfi_bounds(N) ./ N),
                color     = :black,
                linestyle = :dash,
                linewidth = 1.5)
    end

    Colorbar(fig[1, 2], hm, label=L"$$QFI/N")

    if out in (:d, :sd, :ds)
        return T_values, collect(theta_values), qfi_matrix, fig
    elseif out in (:s, :sd, :ds)
        return save(filename, fig, px_per_unit=4)
    else
        return fig
    end
end

function animate_heatmap_Tr_cpp(n=2, hn=3, thn=3;
                                  exe          = "../build/trun",
                                  filename::String = "results/QFI_hm_N$(3n)_animh.mkv",
                                  verbose::Bool = true)

    N = 3n
    h_values= range(0.0, 2.0, length=2^hn+1)
    theta_values = range(0, pi/2, length=2^thn+1)
    T_values     = exp.(range(log(1e-4), log(3.0), length=2^(thn+2)+1))

    tick_ratios  = [0, 0.5, 1, 2, Inf]
    tick_pos     = atan.(tick_ratios)
    tick_labels  = ["0", "0.5", "1", "2", "∞"]
    log_T_values = log10.(T_values)
    log_T_ticks  = floor(log_T_values[1]):1:ceil(log_T_values[end])

    mkpath("results")

    #------------------------------------------------------------------
    # Bootstrap: get qfi_matrix for first h to build the static figure
    #------------------------------------------------------------------
    verbose && println("Creating first frame (h = $(round(h_values[1], digits=2)))...")
    _, _, qfi_matrix_init, _ = hm_Tr_cpp(n, thn, h_values[1];
                                               exe          = exe,
                                               out          = :d,
                                               verbose      = false)

    #------------------------------------------------------------------
    # Build the static figure skeleton
    #------------------------------------------------------------------
    fig = Figure(size=(500, 400))
    ax  = Axis(fig[1, 1],
               xlabel = L"J_2/J_1",
               ylabel = L"$T$ (in log scale)",
               xticks = (tick_pos, tick_labels),
               yticks = (log_T_ticks, [L"10^{%$(Int(e))}" for e in log_T_ticks]))
    xlims!(ax, (theta_values[1], theta_values[end]))
    ylims!(ax, (log_T_values[1], log_T_values[end]))

    mat_obs   = Observable(Float64.(qfi_matrix_init))
    htext_obs = Observable(@sprintf("h=%.2f", h_values[1]))

    hm = heatmap!(ax,
                  collect(theta_values),
                  log_T_values,
                  mat_obs,
                  colormap   = COLOUR,
                  colorrange = (0.0, 5.0))

    text!(ax, 0.02, 0.98,
          text  = htext_obs,
          space = :relative,
          color = :white,
          align = (:left, :top))

    contour!(ax,
             collect(theta_values),
             log_T_values,
             mat_obs,
             levels    = Float64.(qfi_bounds(N) ./ N),
             color     = :black,
             linestyle = :dash,
             linewidth = 1)

    Colorbar(fig[1, 2], hm, label=L"$$QFI/N")

    #------------------------------------------------------------------
    # Record — one frame per h value
    #------------------------------------------------------------------
    p = verbose ? Progress(length(h_values); dt=0.5, barlen=50, color=:yellow) : nothing

    record(fig, filename, h_values; framerate=1) do h

        _, _, qfi_matrix, _ = hm_Tr_cpp(n, thn, Float64(h);
                                              exe          = exe,
                                              out          = :d,
                                              verbose      = false,
                                              clabels      = false)

        mat_obs[]   = Float64.(qfi_matrix)
        htext_obs[] = @sprintf("h=%.2f", h)

        verbose && next!(p)
    end

    verbose && println("Animation saved to $filename")
    return fig
end
