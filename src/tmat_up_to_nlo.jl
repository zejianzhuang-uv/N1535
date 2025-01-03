

using LinearAlgebra
using NLsolve
using Statistics

include("potential.jl")
include("loop_function.jl")
include("phase_space.jl")


function tdet_up_to_nlo_degen(w, par, C, D, L, a; mu=630e0, rs=:rs21, n=2)
    id = Matrix{Float64}(I, n, n)
    wt = WT(w, par, C, n=n)
    vnlo = vnlo_swave(w, par, D, L, n=n)
    gl = Gdr_mat_degenerate(w, par, a, mu=mu, rs=rs, n=n)
    v = wt + vnlo
    return det(id - v * gl)
end


function plot_tdet_up_to_nlo_degen(ax, rew, imw, par, C, D, L, a; rs=:rs21, mu=630e0, n=2)
    dim = length(rew)
    abstdet = zeros(dim, dim)
    for i in 1:dim
        for j in 1:dim
            w = rew[j] + imw[i] * 1im
            abstdet[i, j] = abs(tdet_up_to_nlo_degen(w, par, C, D, L, a, mu=mu, rs=rs, n=n) )
        end
    end
    ax.set_ylim(minimum(imw), maximum(imw))
    ax.set_xlim(minimum(rew), maximum(rew))
    ax.set_xlabel(L"Re$[E_\mathrm{cm}]$ [MeV]")
    ax.set_ylabel(L"Im$[E_\mathrm{cm}]$ [MeV]")
    l = range(minimum(abstdet), maximum(abstdet), 20)
    #ax.contour(rew, imw, abstdet, colors=["#000", "#000"], linestyles="--", levels=l, linewidths=0.9)
    cf = ax.contourf(rew, imw, abstdet, levels=l, cmap=:turbo)
    ax.grid(false)
    return cf
end


function tmat_up_to_nlo_swave(w, par, C, D, L, a; mu=630e0, n=6, ch=:ch11)
    id = Matrix{Float64}(I, n, n)
    # wt = WT(w, par, C, n=6)
    # vnlo = vnlo_swave(w, par, D, L, n=6)
    v = V_up_to_nlo_swave(w, par, C, D, L, n=6)
    # v = wt + vnlo
    gl = Gdr_mat(w, par, a, mu=mu)
    t = inv(id - v * gl) * v

    if ch == :ch11
        return t[1, 1]
    elseif ch == :ch12
        return t[1, 2]
    elseif ch == :ch13
        return t[1, 3]
    elseif ch == :ch14
        return t[1, 4]
    elseif ch == :ch15
        return t[1, 5]
    elseif ch == :ch16
        return t[1, 6]
    end
end

function dist_up_to_nlo_swave(w, par, C, D, L, a, ch)
    hbarc_sq = 0.389379372e6 # MEV × mb
    t = tmat_up_to_nlo_swave(w, par, C, D, L, a, ch=ch)
    psp = phase_space(w, par, ch)
    return psp * abs2(t) * hbarc_sq
end

