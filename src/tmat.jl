

using LinearAlgebra
using NLsolve
import Statistics

include("potential.jl")
include("loop_function.jl")
include("subtraction_cons.jl")
include("phase_space.jl")



function tmat_lo_degenrate(w, par, C, a; mu=630.0f0, rs=:rs21, n=2, ch=:ch11)
    id = diagm(ones(n))
    wt = WT(w, par, C, n=n)
    gl = Gdr_mat_degenerate(w, par, a, mu=mu, rs=rs, n=n)

    t = inv(id - wt * gl) * wt
    if ch == :ch11
        return t[1, 1]
    elseif ch == :ch22
        return t[2, 2]
    end
end

function tdet_lo_degenrate(w, par, C, a; mu=630.0f0, rs=:rs21, n=2)
    id = diagm(ones(n))
    wt = WT(w, par, C, n=n)
    gl = Gdr_mat_degenerate(w, par, a, mu=mu, rs=rs, n=n)
    return det(id - wt * gl)
end

function pole_lo_degen(init_x::Vector, par, C, a; mu=630.0f0, rs=:rs21, n=2)
    sol = nlsolve(x -> cmplx(tdet_lo_degenrate(x[1] - x[2] * 1im, par, C, a, mu=mu, rs=rs)), init_x)
    #return sol
    if sol.residual_norm <= 1e-7
        x = sol.zero
        return x[1] + x[2] * 1im
    else
        return missing
    end

end


function plot_tdet_lo_degenrate(ax, rew, imw, par, C, a; rs=:rs21, mu=630.0f0)
    dim = length(rew)
    abstdet = zeros(dim, dim)
    for i in 1:dim
        for j in 1:dim
            w = rew[j] + imw[i] * 1im
            abstdet[i, j] = abs(tdet_lo_degenrate(w, par, C, a, mu=mu, rs=rs) )
        end
    end
    ax.set_ylim(minimum(imw), maximum(imw) )
    ax.set_xlim(minimum(rew), maximum(rew) )
    ax.set_xlabel(L"Re$[E_\mathrm{cm}]$ [MeV]")
    ax.set_ylabel(L"Im$[E_\mathrm{cm}]$ [MeV]")
    l = range(minimum(abstdet), maximum(abstdet), 20)
    ax.contour(rew, imw, abstdet, colors=["#000", "#000"], linestyles="--", levels=l, linewidths=0.9)
    cf = ax.contourf(rew, imw, abstdet, levels=l, cmap=:turbo)
    ax.grid(false)
    return cf
end

function tmat_lo(w, par, C, a; mu=630e0, n=6, ch=:ch11, born=true)
    id = diagm(ones(n))
    wt = V_lo(w, par, C, n=6, born=born)
    gl = Gdr_mat(w, par, a, mu=mu)
    t = inv(id - wt * gl) * wt
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




function dist_lo(w, par, C, a, ch; born=true)
    hbarc_sq = 0.389379372e6 # MEV × mb
    t = tmat_lo(w, par, C, a, ch=ch, born=born)
    psp = phase_space(w, par, ch)
    return psp * abs2(t) * hbarc_sq
end



function dist_lo_error_band(w, par, C, qmax_sample, ch; mu=630e0)
    dist_w = []
    for ww in w
        dist_a = []
        for qmax in qmax_sample
            a = [subtraction_cons(mchi, mu, qmax) for mchi in par[:mch]]
            dist = dist_lo(ww, par, C, a, ch)
            push!(dist_a, dist)
        end
        push!(dist_w, dist_a)
    end
    return sqrt.( diag(Statistics.cov(hcat(dist_w...)) ) )
end
