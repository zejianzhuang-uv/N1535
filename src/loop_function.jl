

using LinearAlgebra



include("formula.jl")


function Gdr_mat(w, par::AbstractDict, a; mu::Float64=630e0)
    mch = par[:mch]
    gl = diagm([Gdr(w, mchi, aa, mu=mu) for (mchi, aa) in zip(mch, a) ] )
    return gl
end


function Gdr_mat_degenerate(w, par, a::Tuple; mu=630e0, rs=:rs21, n=2)
    mch = par[:mch][1:n]
    if rs == :rs11 && n == 2
        gl = diagm([Gdr(w, mchi, aa, mu=mu) for (mchi, aa) in zip(mch, a[1:n])])
    elseif rs == :rs21 && n == 2
        RS = (:rs2, :rs1)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    elseif rs == :rs12 && n == 2
        RS = (:rs1, :rs2)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    elseif rs == :rs22 && n == 2
        RS = (:rs2, :rs2)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    elseif rs == :rs2111 && n == 4
        RS = (:rs2, :rs1, :rs1, :rs1)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    elseif rs == :rs1111 && n == 4
        RS = (:rs1, :rs1, :rs1, :rs1)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    elseif rs == :rs2221 && n == 4
        RS = (:rs2, :rs2, :rs2, :rs1)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    elseif rs == :rs2222 && n == 4
        RS = (:rs2, :rs2, :rs2, :rs2)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    elseif rs == :rs2211 && n == 4
        RS = (:rs2, :rs2, :rs1, :rs1)
        gl = diagm([Gdr(w, mchi, aa, mu=mu, rs=rss) for (rss, mchi, aa) in zip(RS, mch, a[1:n])])
    else
        println("Please check your RS and n setting.")
    end
    return gl
end










"""DR scheme"""
function Gdr(w, mchi, a; mu=630f0, rs::Symbol=:rs1)
    m, M = mchi
    d = M^2 - m^2
    s = w^2
    q = qcm(w, mchi...)

    g = 2 / (16 * π^2) * M * (a + 2log(M / mu) + (-d + s) / (2s) * 2log(m / M) + q / w * (log(s - d + 2q * w) + log(s + d + 2q * w) - log(-s + d + 2q * w) - log(-s - d + 2q * w)))

    if rs == :rs1
        return g
    elseif rs == :rs2
        s = 2im * M * qcm(w, mchi...) / ((4π) * w)
        return g + s
    else
        println("rs should be :rs1 or :rs2.")
    end
end






