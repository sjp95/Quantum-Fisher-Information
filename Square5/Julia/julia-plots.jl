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

function plot_QFI_vs_T(n=2, thn=7, theta::Float64=pi/4, h::Float64=0.0;
                        exe      = "../build/mt",
                        verbose::Bool = true)

    N    = 3n
    J1   = cos(theta)
    J2   = sin(theta)
    hz   = h

    T_values = exp.(range(log(1e-4), log(3.0), length=2^(thn+2)+1))
    Tmin = T_values[1]
    Tmax = T_values[end]
    nT   = length(T_values)

    fmt(x)     = lpad(round(Int, x * 1e6), 8, '0')
    fmt_arg(x) = round(x, digits=6)
    rnd(x)     = round(x, digits=2)

    #------------------------------------------------------------------
    # HDF5 cache
    #------------------------------------------------------------------
    cachefile = "results/cache_QFI_vs_T_N$(N)_th$(fmt(theta))_h$(fmt(h))_nT$(nT).h5"
    mkpath("results")

    qfi_vals = nothing

    if isfile(cachefile)
        try
            verbose && println("Cache found — loading from $cachefile")
            h5open(cachefile, "r") do f
                qfi_vals = read(f, "qfi_values")
                T_cached = read(f, "T_values")
                if length(qfi_vals) == nT
                    verbose && println("Cache loaded successfully.")
                    T_values = T_cached
                else
                    @warn "Cache dimension mismatch — recomputing."
                    qfi_vals = nothing
                end
            end
        catch e
            @warn "Cache read failed: $e — recomputing."
            qfi_vals = nothing
        end
    end

    #------------------------------------------------------------------
    # Cache miss — call C++
    #------------------------------------------------------------------
    if isnothing(qfi_vals)
        cmd = `$exe $N $(fmt_arg(J1)) $(fmt_arg(J2)) $(fmt_arg(hz)) $(fmt_arg(Tmin)) $(fmt_arg(Tmax)) $nT 3`
        verbose && println("Running: $cmd")
        try
            run(pipeline(cmd, stdout=devnull, stderr=devnull))
        catch e
            error("C++ call failed: $e")
        end

        datfile = "../Data/QFI/QFI_$(N)_$(fmt(J1))_$(fmt(J2))_$(fmt(hz)).dat"
        isfile(datfile) || error("Missing output file: $datfile")

        lines = readlines(datfile)
        isempty(lines) && error("Empty output file: $datfile")

        parsed = Dict{Float64,Float64}()
        for l in lines
            isempty(strip(l)) && continue
            cols = parse.(Float64, split(strip(l)))
            parsed[cols[1]] = cols[5]
        end

        cpp_T_vals = collect(keys(parsed))
        qfi_vals   = [parsed[cpp_T_vals[argmin(abs.(cpp_T_vals .- Tv))]] for Tv in T_values]

        try
            h5open(cachefile, "w") do f
                write(f, "qfi_values", qfi_vals)
                write(f, "T_values",   collect(T_values))
            end
            verbose && println("Cache written to $cachefile")
        catch e
            @warn "Cache write failed: $e"
        end
    end

    #------------------------------------------------------------------
    # Plot
    #------------------------------------------------------------------
    legalpha = 0.9
    fig = Figure(size=(350, 250))
    ax  = Axis(fig[1, 1],
               xlabel = L"$T$",
               ylabel = L"QFI/$N$",
               xscale = log10,
#                xticksmirrored=true,
#                yticksmirrored=true,
               xgridcolor  = (:grey, 0.1),
               ygridcolor  = (:grey, 0.1))
#     ax2 = Axis(fig[1,2])
    lines!(ax, T_values, qfi_vals, color=CMAP[0.6], linewidth=1.5)

    # Draw QFI bound lines
    for (i, bound) in enumerate(reverse(qfi_bounds(N)) ./ N)
        bound == 0.0 && continue
        rem(i,2) == 0 && continue
        hlines!(ax, bound,
                linestyle = LSTYLES[rem(length(qfi_bounds(N))-i+1,length(LSTYLES)-1)+2],
                color     = :black,
                linewidth = 1,
                alpha     = 0.5,
                label     = L"$k=%$(length(qfi_bounds(N))-i)$",
                )
    end
#     for (i, bound) in enumerate(qfi_bounds(N) ./ N)
#         bound == 0.0 && continue
#         rem(i,2) == 0 && continue
#         hlines!(ax, bound,
#                 linestyle = LSTYLES[rem(i,length(LSTYLES)-1)+2],
#                 color     = :black, alpha=0.6,
#                 linewidth = 1)
#     current_limits = ax.finallimits[]
#     ymax = current_limits.origin[1] + current_limits.widths[1]
#     println(ymax)
#     i%2==0 && (text!(ax, 4, bound,
#           text  = L"k=%$(i-1)",
# #           space = :relative,          # x is in relative coords (0–1), y in data coords
#           align = (:left, :center),
#           fontsize = 11,
#           clip_planes = Makie.Plane3f[],
#           color = :black))
#     end
    fig[1,2] = Legend(fig, ax, "Bound for", labelsize=12,
#         backgroundcolor=(:white,legalpha),
        framevisible=false)
#     axislegend(ax2, L"$$Bound for", position=:rt, labelsize=12, backgroundcolor=(:white,legalpha))

    filename = "results/QFI_vs_T_N$(N)_th$(rnd(theta))_h$(rnd(h)).png"
    save(filename, fig, px_per_unit=4)
    verbose && println("Saved to $filename")

    return fig
end

function hm_Tr_cpp(n=2, thn=3, h::Float64=0.0;
                        exe          = "../build/mt",
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
        cmd = `$exe $N $(fmt_arg(J1)) $(fmt_arg(J2)) $(fmt_arg(hz)) $(fmt_arg(Tmin)) $(fmt_arg(Tmax)) $nT 3`
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

    fig = Figure(size=(350, 250))
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
                                  exe          = "../build/mt",
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

function hm_h_cpp(n=2, thn=3, hn=3;
                  exe          = "../build/st",
                  out::Symbol  = :f,
                  verbose::Bool = true,
                  clabels::Bool = true)

    N            = 3n
    theta_values = range(0, pi/2, length=2^thn + 1)
    h_values     = range(0.0, 2.0, length=2^hn  + 1)
    filename     = "results/trimerMag_hm_h_N$(N)_thn$(thn)_hn$(hn).png"

    fmt(x)     = lpad(round(Int, x * 1e6), 8, '0')
    fmt_arg(x) = round(x, digits=6)

    #------------------------------------------------------------------
    # HDF5 cache
    #------------------------------------------------------------------
    cachefile = "results/cache_hm_h_N$(N)_thn$(thn)_hn$(hn).h5"
    mkpath("results")

    sz_matrix = nothing

    if isfile(cachefile)
        try
            verbose && println("Cache found — loading from $cachefile")
            h5open(cachefile, "r") do f
                cached = read(f, "sz_matrix")
                th_c   = read(f, "theta_values")
                h_c    = read(f, "h_values")
                if size(cached) == (length(theta_values), length(h_values))
                    verbose && println("Cache loaded successfully — skipping C++ calls.")
                    sz_matrix      = cached
                    theta_values   = th_c
                    h_values       = h_c
                else
                    @warn "Cache dimensions mismatch — recomputing."
                end
            end
        catch e
            @warn "Cache read failed: $e — recomputing."
        end
    end

    #------------------------------------------------------------------
    # Cache miss — run C++ for each (theta, h) point
    # sz_flag = 4, no qfi, no clean: code = 4
    #------------------------------------------------------------------
    if isnothing(sz_matrix)
        sz_matrix = zeros(length(theta_values), length(h_values))

        # dummy T args — Sz() ignores them but main() requires them
        Tmin = 1e-4
        Tmax = 1e-4
        nT   = 1

        p = verbose ? Progress(length(theta_values) * length(h_values);
                               dt=0.5, barlen=50, color=:yellow) : nothing

        for (ir, theta) in enumerate(theta_values)
            J1 = cos(theta)
            J2 = sin(theta)

            for (ih, h) in enumerate(h_values)
                cmd = `$exe $N 0.001 $(fmt_arg(J1)) $(fmt_arg(J2)) $(fmt_arg(h)) 5`

                try
                    run(pipeline(cmd, stdout=devnull, stderr=devnull))
                catch e
                    @warn "C++ call failed (theta=$(round(theta,digits=3)), h=$(round(h,digits=3))): $e"
                    verbose && next!(p)
                    continue
                end

                datfile = "../Data/Magnetization/Sz/Sz_$(N)_$(fmt(J1))_$(fmt(J2))_$(fmt(h)).dat"

                if !isfile(datfile)
                    @warn "Missing: $datfile"
                    verbose && next!(p)
                    continue
                end

                lines = filter(!isempty ∘ strip, readlines(datfile))
                if isempty(lines)
                    @warn "Empty: $datfile"
                    verbose && next!(p)
                    continue
                end

                # last line — in case file has multiple entries
                cols = parse.(Float64, split(strip(lines[end])))
                sz_matrix[ir, ih] = cols[4]

                verbose && next!(p)
            end
        end

        #--------------------------------------------------------------
        # Write cache
        #--------------------------------------------------------------
        try
            h5open(cachefile, "w") do f
                write(f, "sz_matrix",    sz_matrix)
                write(f, "theta_values", collect(theta_values))
                write(f, "h_values",     collect(h_values))
            end
            verbose && println("Cache written to $cachefile")
        catch e
            @warn "Cache write failed: $e"
        end
    end

    #------------------------------------------------------------------
    # Plot
    #------------------------------------------------------------------
    tick_ratios = [0, 0.5, 1, 2, Inf]
    tick_pos    = atan.(tick_ratios)
    tick_labels = ["0", "0.5", "1", "2", "∞"]

    fig = Figure(size=(400, 350))
    ax  = Axis(fig[1, 1],
               xlabel = L"J_2/J_1",
               ylabel = L"h",
               xticks = (tick_pos, tick_labels),
               yticks = (range(0, 2, step=0.5), ["0", "0.5", "1", "1.5", "2"]))
    xlims!(ax, (theta_values[1], theta_values[end]))
    ylims!(ax, (h_values[1], h_values[end]))

    hm = heatmap!(ax,
                  collect(theta_values),
                  collect(h_values),
                  sz_matrix,
                  colormap   = COLOUR,
                  colorrange = (0.0, 0.5))   # Sz/N in [0, 0.5] for spin-1/2

    Colorbar(fig[1, 2], hm, label=L"\langle S_z \rangle / N")

    if out == :f
        return fig
    elseif out == :s
        save(filename, fig, px_per_unit=4)
        return fig
    elseif out == :d
        return collect(theta_values), collect(h_values), sz_matrix, fig
    end
end

# for h in [0.0 0.2 0.25 0.5 0.7 0.9 1.0 1.5]
#     hm_Tr_cpp(4,7,h; out=:s,clabels=true)
# end
